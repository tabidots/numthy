(ns numthy.helpers
  (:require [clojure.string :as s]
            [clojure.math.numeric-tower :as tower]))

;; GENERAL UTILITIES

(defn pfilter
  "Like filter, except uses multiple cores. Useful for huge collections."
  [pred coll]
  (->> (pmap #(when (pred %) %) coll)
       (remove nil?)))

;; PRIMALITY HELPERS

(defn isqrt
  "floor(√n). When incremented, provides an upper bound for factorization."
  ;; Java interop is super fast but not accurate for n > 1E24 (approx) due to
  ;; floating-point rounding. Uses a slightly slower but pinpoint-precise method for n > 1E24.
  [n]
  (if (< n 1E24)
    (-> (Math/sqrt n) bigint)
    ;; https://cs.stackexchange.com/a/30383
    (let [half-bit-length (quot (.bitLength (bigint n)) 2)]
      (loop [a (tower/expt 2 half-bit-length)
             b a
             c (*' a a)]
        (cond
          (zero? b) a
          (> c n)   (recur (-' a b) (quot b 2) (+ c (*' -2 a b) (*' b b)))
          :else     (recur (+' a b) (quot b 2) (+ c (*' 2 a b) (*' b b))))))))

(defn divisible? [dividend divisor]
  (zero? (mod dividend divisor)))

(defn naive-prime? [n]
  (cond
    (<= n 1)  false
    (= n 2)   true
    (even? n) false ;; Will also weed out non-integers
    :else     (let [lim (inc (isqrt n))]
                (loop [i 3]
                  (cond
                    (>= i lim)       true
                    (divisible? n i) false
                    :else            (recur (+ i 2)))))))

(def prime? (memoize naive-prime?))

;; FACTORIZATION

(defn factors
  [n]
  (->> (range 1 (inc (isqrt n)))
       (filter #(divisible? n %))
       (mapcat (fn [x] [x (/ n x)]))
       (into (sorted-set))))

(defn phi
  "Euler's totient function, optimized using the product rule."
  [n]
  ;; https://en.wikipedia.org/wiki/Euler%27s_totient_function#Example
  (let [prime-factors (filter prime? (factors n))
        diffs         (map #(- 1 (/ 1 %)) prime-factors)]
     (reduce * n diffs)))

(defn prime-factorization [n]
  "Returns the prime factorization of an integer, e.g., 168 -> [2 2 2 3 7]."
  (when (and (integer? n) (pos? n))
    (loop [n       n
           primes  (if (> n 1000000)
                     (filter prime? (factors n))
                     (filter prime? (range)))
           factors []]
      (let [k (first primes)]
        (cond
          (= 1 n)          factors
          (divisible? n k) (recur (/ n k) primes (conj factors k))
          :else            (recur n (rest primes) factors))))))

;; TODO: Reimplement this using hash-maps

(defn prime-powers-integer [n]
  "Returns the exponent vector for the prime power representation of an integer,
  e.g., 168 = 2*2*2*3*7 = 2^3 * 3^1 * 5^0 * 7^1 -> (3 1 0 1)"
  (when (and (integer? n) (pos? n))
    (if (= n 1) [0]
      (let [pf      (prime-factorization n)
            primes  (filter prime? (range (inc (peek pf))))
            freqs   (frequencies pf)
            get-exp (fn [p] (if-let [exp (get freqs p)]
                              exp 0))]
        (map get-exp primes)))))

(defn prime-powers-rational [q]
  (let [r (rationalize q)]
    ;; Sanity check to accept decimal representations of rational numbers
    ;; while still rejecting irrational numbers
    (when (or (ratio? q) (= (double r) (double q)))
      (let [n (prime-powers-integer (numerator r))
            d (prime-powers-integer (denominator r))
            pn (concat n (repeat (count d) 0))
            pd (concat d (repeat (count n) 0))]
        (->> (map - pn pd)
             reverse
             (drop-while zero?) ;; truncate zeros from the end
             reverse)))))

(defn prime-powers [n]
  (when (pos? n)
    (if (integer? n)
      (prime-powers-integer n)
      (prime-powers-rational n))))

(defn prime-powers->num [pp]
  (let [primes (filter prime? (range))]
    (reduce *' (map tower/expt primes pp))))

;; https://en.wikipedia.org/wiki/Prime_omega_function
;; ω(n) = (count (distinct (prime-factorization n)))
;; Ω(n) = (count (prime-factorization n))

(defn smooth? [B n]
  "A positive integer is called B-smooth if none of its prime factors is greater than B. "
  (<= (apply max (prime-factorization n)) B))

(defn regular? [n]
  "5-smooth numbers are also called regular numbers."
  (smooth? 5 n))

(defn humble? [n]
  "7-smooth numbers are also called humble numbers."
  (smooth? 7 n))

(defn powersmooth? [B n]
  "A positive integer is called B-powersmooth if none of its prime powers are greater than B."
  (let [primes (filter prime? (range))]
    (every? #(<= % B) (map tower/expt primes (prime-powers n)))))

(defn radical [n]
  "The product of the distinct prime numbers dividing a positive integer n.
  This is also the largest square-free factor of n."
  (reduce * (distinct (prime-factorization n))))

(defn square-free? [n]
  (or (= 1 n)
      (prime? n) ; If n is a known prime, this is much faster
      (apply distinct? (prime-factorization n))))

(defn powerful? [n]
  "A powerful number is a positive integer n such that for every prime factor p of n,
  p2 is also a factor of n."
  (let [pf (filter prime? (factors n))]
    (every? #(divisible? n (* % %)) pf)))

(defn liouville [n]
  ; https://en.wikipedia.org/wiki/Liouville_function
  (tower/expt -1 (count (prime-factorization n))))

(defn mobius [n]
  ; https://en.wikipedia.org/wiki/M%C3%B6bius_function
  (if (square-free? n)
    (if (even? (count (prime-factorization n))) 1 -1)
    0))

; https://en.wikipedia.org/wiki/Exponentiation_by_squaring
(defn exp-by-sq [x n]
  (cond
    (neg? n)   (exp-by-sq (/ 1 x) (- n))
    (zero? n)  1
    (= 1 n)    x
    (even? n)  (exp-by-sq (* x x) (/ n 2))
    (odd? n)   (exp-by-sq (* x x) (-> (- n 1) (/ 2)))))

(defn factorial [n]
  (reduce * (take n (iterate (partial inc) (bigint 1)))))

(defn chebyshev [x]
  ;; Second Chebyshev ψ(x)
  ;; https://en.wikipedia.org/wiki/Chebyshev_function
  (->> (filter prime? (range (inc x)))
       (map #(*' (Math/floor (/ (Math/log x)
                                (Math/log %)))
                 (Math/log10 %)))
       (reduce +')))

(defn von-mangoldt [n]
  ;; Λ(n)
  ;; https://en.wikipedia.org/wiki/Von_Mangoldt_function
  (let [pf (prime-factorization n)]
    (if (= 1 (count (distinct pf)))
      (Math/log10 (first pf))
      0)))

;; SEQUENCES

(defn fibonacci
  ([]
   (->> (iterate (fn [[f0 f1]] [f1 (+ f0 f1)]) [0N 1N])
        (map peek)))
  ([n]
   (take n (fibonacci))))
