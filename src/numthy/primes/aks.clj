(ns numthy.primes.aks
  (:require [clojure.math.numeric-tower :refer [sqrt expt]]
            [numthy.modular-arithmetic.multiplicative-order :refer [multiplicative-order]]
            [numthy.helpers :refer [divisible?]]
            [numthy.polynomials.core :as p]
            [numthy.perfect-powers :refer [perfect-power?]]
            [numthy.factorization.core :refer [phi]]))

(defn p-every?
  "Multi-core version of `every?`."
  [pred coll]
  (every? true? (pmap pred coll)))

(defn aks-prime?
  "Uses the Agrawal–Kayal–Saxena primality test to determine if an integer n
  is prime. Returns true if prime, nil or false otherwise."
  [n]
  ;; 1. Check if n is a perfect power
  (when-not (or (perfect-power? n) (even? n))
    ;; 2. Find the smallest r such that ord_r(n) > (log_2 n)^2.
    (let [log (-> (Math/log n) (/ (Math/log 2)) (expt 2))
          r   (some (fn [r] (when-let [ord (multiplicative-order n r)]
                              (when (> ord log) r)))
                    (drop 2 (range)))
          lim (min r (dec n))]
      ;; 3. For all 2 ≤ a ≤ min(r, n−1), check that a does not divide n
      ;; (composite if so) [already took care of case where a=2 by returning nil for even n]
      (when (not-any? #(divisible? n %) (range 3 (inc lim) 2))
        ;; 4. If n ≤ r, output prime.
        (if (<= n r) true
          (let [log2n (/ (Math/log n) (Math/log 2))
                lim   (->> (sqrt (phi r)) (* log2n) bigint)
                lhs   (fn [a] (p/poly-rem (p/mod-exp {1 1, 0 a} n n)
                                          {r 1, 0 -1}))
                rhs   (fn [a] (p/poly-rem {n 1, 0 a}
                                          {r 1, 0 -1}))]
            ;; 5. If (X+a)^n != (X^n)+a (mod X^r − 1,n) for ANY a from 1 to lim,
            ;; n is composite.
            ;; In other words, prime? = true iff (X+a)^n = (X^n)+a (mod X^r − 1,n)
            ;; for ALL a from 1 to lim
            (p-every? (fn [a] (= (lhs a) (rhs a))) (range 1 lim))))))))

;; TODO: Final step of AKS algorithm is nearly intractable, even for numbers that
;; aren't that large. How to optimize it?
