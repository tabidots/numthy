(ns numthy.modular-arithmetic.multiplicative-order
  (:require [clojure.math.numeric-tower :refer [expt gcd lcm]]
            [numthy.factorization.core :refer [factorize prime-factors-with-multiplicity]]
            [numthy.helpers :refer [coprime?]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]
            [numthy.modular-arithmetic.primitive-roots :refer [powers-of-a-mod-n]]))

(comment
  (defn naive-multiplicative-order
    "Find the smallest positive integer k such that a^k ≡ 1 (mod n) where
    a is coprime to n. a^0 ≡ 1 (mod n) for any n. "
    [a n]
    (when (coprime? a n)
      (first (filter #(= 1 (mod-pow a % n)) (rest (range))))))

  "For multiplicative order, an intuitive way to find k for small n is to count
  the distinct powers of a (mod n). However, this is clearly intractable for most n."
  (defn naive-multiplicative-order-2
    [a n]
    (when (coprime? a n)
      (count (distinct (powers-of-a-mod-n a n))))))

;; Actual implementation

(defn- ord-prime
  "If a is not a primitive root of a prime p, then ord_p(a) is the smallest remainder
  of φ(p)/x, where x is a divisor of φ(p), s.t. a^(φ(p)/x) ≡ 1 mod p.
  If it is a primitive root, then ord_p(a) = φ(p)."
  [a p]
  (cond
    (= (mod a p) 1)       1
    (= (mod a p) (dec p)) 2 ;; a ≡ -1 mod p → ord_p(a) = 2
    :else
    (let [phi (dec p)]
      (->> (prime-factors-with-multiplicity phi)
           (reductions / phi)
           (filter #(= 1 (mod-pow a % p)))
           (apply min)))))

(defn- ord-prime-power
  "Finds ord_q(a), where q is a prime power of the form p^n."
  ;; Translated from orderpn(a,p,n) in this http://www.numbertheory.org/gnubc/phi
  [a p n]
  (when (and (> a 1) (coprime? a p))
    (if (= n 1) (ord-prime a p) ;; Sanity check
      (let [d (if (= p 2) 2 (ord-prime a p))
            h (->> (rest (range (inc n)))
                   (take-while #(= 1 (mod-pow a d (expt p %))))
                   last)] ;; The largest 1 ≦ h ≦ n s.t. a^d ≡ 1 mod p^h
        (if (= n h) d
          (* d (expt p (- n h))))))))

(defn multiplicative-order
  "Finds ord_m(a), where a > 1 and coprime to m, which is the smallest positive
  integer k such that a^k ≡ 1 (mod n) where a is coprime to n. Returns nil for
  a = 0 and a = 1."
  ;; Translated from orderm(a,m) in this http://www.numbertheory.org/gnubc/phi
  ;; By the Chinese Remainder Theorem, it's enough to calculate the lcm of the
  ;; multiplicative orders of a for each prime exponent p^k of m. (From Rosetta Code)
  [a m]
  (when (and (> a 1) (coprime? a m))
    (reduce-kv (fn [o p n]
                 (lcm o (ord-prime-power a p n)))
               1 (factorize m))))
