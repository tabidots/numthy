(ns numthy.modular-arithmetic.primitive-roots
  (:require [clojure.math.numeric-tower :refer [expt]]
            [flatland.useful.seq :refer [take-until]]
            [numthy.factorization.core :refer [distinct-prime-factors divide-out-small-factors phi]]
            [numthy.helpers :refer [coprime? rand-num]]
            [numthy.modular-arithmetic.groups :refer [multiplicative-group cyclic?]]
            [numthy.modular-arithmetic.utils :refer [mod-mul mod-pow odd-prime? odd-prime-power?]]))

(defn- powers-of-a-mod-n*
  "a^k (mod n) for all 0 ≤ k < n, where k ∈ ℤ."
  [a n]
  (map #(mod-pow a % n) (range n)))

(defn- take-while-distinct [coll]
  ;; https://github.com/taylorwood/advent-of-code/blob/master/src/advent_of_code/2017/16.clj
  (let [seen (volatile! #{})]
    (take-while (fn [v]
                  (if (@seen v) false
                    (vswap! seen conj v)))
                coll)))

(defn powers-of-a-mod-n
  "a^k (mod n), where k ∈ ℤ, until a duplicate is found. For prime moduli, the effective
  range is still 0 ≤ k < n. However, for composite moduli, this will short-circuit at the
  beginning of the second cycle. Note that the duplicate is included in the sequence."
  [a n]
  (if (odd-prime? n)
    (powers-of-a-mod-n* a n)
    (take-while-distinct (powers-of-a-mod-n* a n))))

;; TODO: Use batch GCD here?

(defn naive-primitive-roots
  "Primitive roots of (ℤ/nℤ,×), the multiplicative group of integers mod n,
  i.e., all of the integers < n that, when multiplied repeatedly by themselves,
  generate the entire group. This is essentially a brute-force method
  that will take a long time for n > 500 or so."
  [n]
  (let [group (set (multiplicative-group n))
        roots (filter #(= group (set (powers-of-a-mod-n % n))) group)]
    (when (not-empty roots)
      (into (sorted-set) roots))))

(defn primitive-root?
  "Efficiently determines if x is a primitive root mod n by leveraging the fact
  that ord_n(x) will always be Φ(n) or one of its divisors. If there is no prime
  divisor of Φ(n) s.t. x^(Φ(n)/d) ≡ 1 mod n, then ord_n(x) must be Φ(n) mod n,
  which means that x is a primitive root mod n."
  ; https://stackoverflow.com/a/26636457/4210855
  ([x n]
   (when (< 1 x n)
     (when-let [phi' (cond
                       (odd-prime? n)       (dec n)
                       (and (cyclic? n)
                            (coprime? x n)) (phi n)
                       :else                nil)]
       (let [pfs (distinct-prime-factors phi')]
         (primitive-root? x n phi' pfs)))))
  ([x n phi' pfs]
   (not-any? #(= 1 (mod-pow x (/ phi' %) n)) pfs)))

(defn least-primitive-root
  "Given a prime p, finds the least integer 2 ≤ a < p, s.t. a^φ(p) ≡ 1 mod p."
  ;; https://math.stackexchange.com/a/133720
  ;; But don't waste time factoring φ(p) completely, because the LPR is not likely to be very large
  [p]
  (if (= p 2) 1
    (when (odd-prime? p)
      (let [phi' (phi p)
            pfs (distinct (keys (:factors (divide-out-small-factors phi'))))]
        (->> (drop 2 (range))
             (filter #(primitive-root? % p phi' pfs))
             first)))))

(defn random-primitive-root
  "Given a prime p > 2, generates random numbers until it finds a primitive rood mod p."
  [p]
  (when (odd-prime? p)
    (let [phi' (phi p) pfs (distinct-prime-factors phi')]
      (->> #(rand-num 2 p)
           repeatedly
           (filter #(primitive-root? % p phi' pfs))
           first))))

(defn primitive-root-probability
  "The probability of choosing a primitive root mod n (n > 2) at random."
  [n]
  (when (> n 2)
    (let [phi' (phi n)]
      (/ (phi phi') phi'))))

(defn primitive-roots-prime
  "Uses the least primitive root of an odd prime p to find all other primitive roots."
  ;; https://math.stackexchange.com/a/133720
  ;; if a is a primitive root mod p, a can generate all other remainders 1…(p−1) as powers.
  ;; a^m mod p is another primitive root iff m and p−1 are coprime.
  [p]
  (let [lpr (least-primitive-root p) phi' (phi p)]
    (into (sorted-set)
          (comp (filter #(coprime? % phi'))
                (map #(mod-pow lpr % p)))
          (range 1 phi'))))

(defn- primitive-roots-prime-squared
  ;; http://courses.mai.liu.se/GU/TATA54/Lectures/lecture5.pdf
  [p]
  (let [n    (*' p p)
        lift (fn [a] (map #(+ a (* % p)) ;; a+tp for 0 ≤ t < p
                          (range p)))]
    ;; if a is a p.r. mod p, a+t*p is a p.r. mod p^2 for all 0 ≤ t < p EXCEPT one
    (reduce (fn [res r]
              (->> (lift r)
                   (remove #(= 1 (mod-pow % (dec p) n))) ;; Remove the exception
                   (into res)))
            (sorted-set) (primitive-roots-prime p))))

(defn primitive-roots-prime-power
  ;; https://exploringnumbertheory.wordpress.com/2013/11/01/primitive-roots-of-powers-of-odd-primes/
  ;; http://courses.mai.liu.se/GU/TATA54/Lectures/lecture5.pdf
  [[p k]]
  (cond
    (= k 1)                 (primitive-roots-prime p)
    (and (= p 2) (> k 2))   nil
    :else
    (letfn [(lift [a e] (map #(+ a (* % (expt p e))) ;; a+t*p^e for 0 ≤ t < p
                             (range p)))]
      (loop [j 2 roots (primitive-roots-prime-squared p)]
        (if (= j k) roots
          (recur (inc j)
                 (reduce (fn [res r]                 ;; if a is a p.r. mod p^k, k ≧ 2
                           (into res (lift r j)))    ;; then a+tp is a p.r. mod p^k+1
                         (sorted-set) roots)))))))   ;; if (a+tp)^(p-1) ≢ 1 mod p^k+1

(defn primitive-roots-2pk
  ;; https://exploringnumbertheory.wordpress.com/2013/11/02/primitive-roots-of-twice-the-powers-of-odd-primes/
  [[p k]]
  (when-let [roots (primitive-roots-prime-power [p k])]
    (let [pk (expt p k) n (* 2 pk)]
      (into (sorted-set)                 ;; if a is a p.r. mod p^k
            (map #(if (odd? %) %         ;; and a is odd, then a is a p.r. mod 2*p^k;
                    (mod (+ % pk) n)))   ;; if a is even, a+p^k mod 2*p^k is a p.r. mod 2*p^k
            roots))))

(defn primitive-roots
  [n]
  (cond
    (< n 1) nil
    (= n 1) #{0}
    (= n 2) #{1}
    (= n 4) #{3}
    :else
    (let [opp       (odd-prime-power? n)
          twice-opp (odd-prime-power? (* n 1/2))]
      (cond
        (some? opp)       (primitive-roots-prime-power opp)
        (some? twice-opp) (primitive-roots-2pk twice-opp)
        :else             nil))))
