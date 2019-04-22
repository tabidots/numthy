(ns numthy.p-adic
  (:require [clojure.math.numeric-tower :refer [expt]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.helpers :refer [divisible?]]))

(defn p-adic-order
  "nu_p(n) is the largest exponent v of a prime p such that p^v divides n."
  [p n]
  (if (zero? n) ##Inf
    (when (prime? p)
      (last (take-while (fn [v] (divisible? n (expt p v)))
                        (rest (range)))))))

;; TODO: Extend into the rationals
;; https://en.wikipedia.org/wiki/P-adic_order#Rational_numbers
