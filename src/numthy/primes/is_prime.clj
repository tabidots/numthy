(ns numthy.primes.is-prime
  (:require [numthy.helpers :refer [isqrt divisible?]]))

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

;; TODO: https://en.wikipedia.org/wiki/Quadratic_Frobenius_test
;; TODO: https://en.wikipedia.org/wiki/Pocklington_primality_test
;; TODO: https://en.wikipedia.org/wiki/Lucas_primality_test
;; TODO: https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test
