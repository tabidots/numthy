(ns numthy.perfect-powers
  (:require [clojure.math.numeric-tower :refer [expt]]
            [numthy.helpers :refer [isqrt]]))

(defn babylonian-root
  "High-precision BigDecimal nth-root using the Babylonian algorithm,
  with a close initial approximation for ridiculously fast convergence."
  [n root]
  (cond
    (neg? n)     nil
    (zero? n)    0
    (= 1 n)      1
    (zero? root) nil
    (= 1 root)   n
    :else
    (let [eps 0.000000000000000000000000000000000000000001M]
      (loop [t (bigdec (expt n (/ root)))] ;; rough initial approx
        (let [ts  (repeat (dec root) t)
              nxt (with-precision 100 (-> (/ n (reduce *' ts))
                                          (+ (reduce +' ts))
                                          (/ root)))]
          (if (< (.abs (- nxt t)) eps) t
            (recur nxt)))))))

(defn nth-root-is-integer?
  "Tests if the nth root of x is an integer in the mathematical
  (not programming) sense—i.e., if it is a whole number."
  [x n]
  (when-some [root (some-> (babylonian-root x n) bigint)]
    (let [exp (bigint (reduce *' (repeat n root)))]
      (when (= x exp) {root n}))))

(defn perfect-square?
  "Tests if x is a perfect square."
  [x]
  (when-some [root (isqrt x)]
    (when (= (*' root root) x)
      {root 2})))

(defn perfect-power?
  "If a^b = n for any pair of integers a, b > 1, returns [a b] for the smallest prime b,
  else nil. Excluding even values of b > 2 is for computational efficiency."
  [n]
  (let [max-power (/ (Math/log n) (Math/log 2)) ; https://cr.yp.to/papers/powers-ams.pdf
        powers    (cons 2 (filter odd? (range 3 (inc max-power))))]
    (some #(nth-root-is-integer? n %) powers)))

;; Why is log2(n) the max power? Because 2 is the smallest possible "a" in n = a^b,
;; and log2(a) is the "b" for that "a".

;; Strictly speaking, "powers" above this should be prime powers only, but
;; that would require additional overhead by filtering primes,
;; which seems circuitous since one purpose of this function is already to serve as
;; the first step in an algorithm to test the primality of a gigantic integer
;; Plus, removing the prime? filter allows this function to
;; work without specifying a custom prime? function.

;; Some test values
;; 6221821273427820544 = 2494357888^2, 484^7
;; 92709463147897837085761925410587 = 3^67

(defn perfect-powers
  "Returns a map of keys a > 1 and vals b > 1 such that n = a^b, where a, b ∈ ℤ.
  Does not exclude even powers. Sorted in ascending order of a. Returs nil if n is not
  a perfect power."
  [n]
  (let [max-power (/ (Math/log n) (Math/log 2)) ; https://cr.yp.to/papers/powers-ams.pdf
        powers    (range 2 (inc max-power))
        p-powers  (keep #(nth-root-is-integer? n %) powers)]
    (when (not-empty p-powers)
      (into (sorted-map) p-powers))))
