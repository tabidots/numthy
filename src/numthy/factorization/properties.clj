(ns numthy.factorization.properties
  (:require [clojure.math.numeric-tower :refer [expt]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.helpers :refer [divisible?]]
            [numthy.factorization.core :refer [distinct-prime-factors factorize prime-factors-with-multiplicity prime-signature]]))

(defn smooth? [B n]
  "A positive integer is called B-smooth if none of its prime factors is greater than B. "
  (<= (apply max (distinct-prime-factors n)) B))

(defn regular? [n]
  "5-smooth numbers are also called regular numbers."
  (smooth? 5 n))

(defn humble? [n]
  "7-smooth numbers are also called humble numbers."
  (smooth? 7 n))

(defn powersmooth? [B n]
  "A positive integer is called B-powersmooth if none of its prime powers are greater than B."
  (every? (fn [[p e]] (<= (expt p e) B))
          (factorize n)))

(defn radical [n]
  "The product of the distinct prime numbers dividing a positive integer n.
  This is also the largest square-free factor of n."
  (reduce *' (distinct-prime-factors n)))

(defn square-free? [n]
  (or (= 1 n)
      (prime? n) ;; If n is a known prime, this is much faster
      (every? #(= % 1) (prime-signature n))))

(defn powerful? [n]
  "A powerful number is a positive integer n such that for every prime factor p of n,
  p^2 is also a factor of n."
  (or (= 1 n)
      (every? #(> % 1) (prime-signature))))

(defn liouville [n]
  ; https://en.wikipedia.org/wiki/Liouville_function
  (expt -1 (count (prime-factors-with-multiplicity))))
