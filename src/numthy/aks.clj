(ns numthy.aks
  (:require [clojure.math.numeric-tower :as tower]
            [numthy.modular-arithmetic :as ma]
            [numthy.helpers :as h]
            [numthy.polynomial :as p]))

(defn babylonian-root
  "High-precision BigDecimal nth-root using the Babylonian algorithm,
  with a close initial approximation for ridiculously fast convergence."
  [n root]
  (if (zero? n) 0
    (let [eps 0.000000000000000000000000000000000000000001M]
      (loop [t (bigdec (tower/expt n (/ root)))] ;; rough initial approx
        (let [ts  (repeat (- root 1) t)
              nxt (with-precision 100 (-> (/ n (reduce *' ts))
                                          (+ (reduce +' ts))
                                          (/ root)))]
          (if (< (.abs (- nxt t)) eps) t
            (recur nxt)))))))

(defn nth-root-is-integer?
  "Tests if the nth root of x is an integer in the mathematical
  (not programming) sense—i.e., if it is a whole number.)."
  [x n]
  (let [floor (bigint (Math/floor (babylonian-root x n)))
        exp   (bigint (reduce *' (repeat n floor)))]
    (= x exp)))

;; Some test values
;; 6221821273427820544 = 484^7
;; 92709463147897837085761925410587 = 3^67

(defn perfect-power? [n]
  (let [max-power (/ (Math/log n) (Math/log 2)) ; https://cr.yp.to/papers/powers-ams.pdf
        powers    (cons 2 (filter odd? (range 3 (inc max-power))))]
    (some (partial nth-root-is-integer? n) powers)))

;; Why is log2(n) the max power? Because 2 is the smallest possible "a" in n = a^b,
;; and log2(a) is the "b" for that "a".

;; Strictly speaking, "powers" above this should be prime powers only, but
;; that would require additional overhead by filtering primes
;; which seems circuitous since the purpose of this function is already
;; the first step in an algorithm to test the primality of a gigantic
;; integer. Plus, removing the prime? filter allows this function to
;; work without specifying a custom prime? function.
;;
;; Otherwise, use (filter prime? (range 2 (inc max-power)))

(defn naive-phi
  "Naive version of Euler's totient function that only uses gcd, since
  the optimized version requires factoring n first."
  [n]
  (count (filter #(= 1 (tower/gcd n %)) (range 1 n))))

(defn aks-prime?
  "Uses the Agrawal–Kayal–Saxena primality test to determine if an integer n
  is prime. Returns true if prime, nil or false otherwise."
  [n]
  ;; 1. Check if n is a perfect power
  (when-not (perfect-power? n)
    ;; 2. Find the smallest r such that ord_r(n) > (log_2 n)^2.
    (let [log (-> (Math/log n) (/ (Math/log 2)) (tower/expt 2))
          r   (first (keep (fn [r] (when-let [ord (ma/multiplicative-order n r)]
                                     (when (> ord log) r)))
                           (range)))
          lim (min r (- n 1))]
      ;; 3. For all 2 ≤ a ≤ min(r, n−1), check that a does not divide n
      ;; (composite if so)
      (when (not-any? (partial h/divisible? n) (range 2 (inc lim)))
        ;; 4. If n ≤ r, output prime.
        (if (<= n r) true
          (let [log2n (/ (Math/log n) (Math/log 2))
                lim   (->> (tower/sqrt (naive-phi r)) (* log2n) bigint)
                lhs   (fn [a] (p/poly-rem (p/mod-exp {1 1, 0 a} n n)
                                          {r 1, 0 -1}))
                rhs   (fn [a] (p/poly-rem {n 1, 0 a}
                                          {r 1, 0 -1}))]
            ;; 5. If (X+a)^n != (X^n)+a (mod X^r − 1,n) for ANY a from 1 to lim,
            ;; n is composite.
            ;; In other words, prime? = true iff (X+a)^n = (X^n)+a (mod X^r − 1,n)
            ;; for ALL a from 1 to lim
            (every? (fn [a] (= (lhs a) (rhs a))) (range 1 lim))))))))
