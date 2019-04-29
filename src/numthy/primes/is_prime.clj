(ns numthy.primes.is-prime
  (:require [numthy.helpers :refer [isqrt divisible?]]))

(defn naive-prime?
  "Simple primality testing by trial division of odd integers ≦ √n."
  [n]
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

(defn quick-prime?
  "Quick probabilistic + deterministic check for use in algorithms where
  many numbers need to be checked rapidly, e.g., factorization algos.
  Java's .isProbablePrime returns bizarre results for small numbers sometimes,
  so this is a sure-fire way to avoid NullPointer exceptions deep in other pipelines."
  [n]
  (when (integer? n)  ;; Sanity check
    (if (< n 1000) (prime? n)
      (.isProbablePrime (biginteger n) 5))))


;; TODO: https://en.wikipedia.org/wiki/Quadratic_Frobenius_test
;; TODO: https://en.wikipedia.org/wiki/Pocklington_primality_test
;; TODO: https://en.wikipedia.org/wiki/Lucas_primality_test
;; TODO: https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test
