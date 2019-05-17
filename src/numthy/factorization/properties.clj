(ns numthy.factorization.properties
  (:require [clojure.math.numeric-tower :refer [expt gcd]]
            [numthy.factorization.core :refer [distinct-prime-factors factorize phi prime-factors-with-multiplicity prime-signature]]
            [numthy.primes.is-prime :refer [quick-prime?]]
            [numthy.helpers :refer [divisible?]]))

(comment
 "The following functions determine some mathematical properties of a number n
 based on its prime factorization. Generally, these properties could be determined
 in the process of factorizing the number, so this is circuitous and inefficient,
 especially for large n. However, the functions are useful to illustrate the meaning
 of these properties in words.")

(defn smooth?
  "A positive integer is called B-smooth if none of its prime factors is greater than B. "
  [B n]
  (<= (apply max (distinct-prime-factors n)) B))

(defn regular?
  "5-smooth numbers are also called regular numbers."
  [n]
  (smooth? 5 n))

(defn humble?
  "7-smooth numbers are also called humble numbers."
  [n]
  (smooth? 7 n))

(defn powersmooth?
  "A positive integer is called B-powersmooth if none of its prime powers are greater than B."
  [B n]
  (every? (fn [[p e]] (<= (expt p e) B))
          (factorize n)))

(defn radical
  "The product of the distinct prime numbers dividing a positive integer n.
  This is also the largest square-free factor of n."
  [n]
  (reduce *' (distinct-prime-factors n)))

(defn square?
  [n]
  (even? (reduce gcd (prime-signature n))))

(defn square-free?
  [n]
  (or (= 1 n)
      (quick-prime? n) ;; If n is a known prime, this is much faster
      (every? #(= % 1) (prime-signature n))))

(defn powerful?
  "A powerful number is a positive integer n such that for every prime factor p of n,
  p^2 is also a factor of n."
  [n]
  (or (= 1 n)
      (every? #(> % 1) (prime-signature n))))

(defn achilles?
  "An Achilles number is one that is powerful but not perfect (a perfect power)."
  [n]
  (when (> n 1)
    (let [sig (prime-signature n)]
      (and (every? #(> % 1) sig) (= 1 (reduce gcd sig))))))

(defn strong-achilles?
  "A strong Achilles number is an Achilles number n for which φ(n) is also an Achilles number."
  [n]
  (and (achilles? n) (achilles? (phi n))))

(defn sphenic?
  "A sphenic number is the product of three distinct primes."
  [n]
  (= (prime-signature n) [1 1 1]))

(defn liouville [n]
  ; https://en.wikipedia.org/wiki/Liouville_function
  (expt -1 (count (prime-factors-with-multiplicity))))
