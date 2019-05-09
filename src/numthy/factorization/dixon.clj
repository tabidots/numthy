(ns numthy.factorization.dixon
  "Implementation of Dixon's integer factorization algorithm, which belongs to the family of congruence-of-squares methods."
  (:require [clojure.core.reducers :as r]
            [numthy.factorization.squares-utils :refer [factorize-rhs find-factor smooth-factorization smoothness-bound]]
            [numthy.factorization.trial-division :refer [batch-smooth-filter]]
            [numthy.helpers :refer [isqrt]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]
            [numthy.primes.sieves :refer [primes-to]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

(defn- add-relations
  "As per Dixon's algorithm, runs up the range [√n, n] testing for values of x
  s.t. x^2 mod n is smooth over the factor base. Smooth relations are added to a
  map of the form {x x^2 mod n}. This function consumes and tests candidates in
  batches of 250."
  [n factor-base limit]
  (->> (partition 250 (range (inc (isqrt n)) (inc n)))
       (r/reduce (fn [smooths xs]
                   (if (> (count smooths) limit) (reduced smooths)
                     (->> (batch-smooth-filter factor-base #(mod-pow % 2 n) xs)
                          (factorize-rhs factor-base)
                          (into smooths))))
                 {})))

(defn dixon
  "Uses Dixon's factorization algorithm to find one prime factor of n.
  Returns nil if algorithm fails to find a factor."
  ;; With reference to https://crypto.stanford.edu/cs359c/17sp/projects/BrendonGo.pdf
  [n]
  (let [B           (smoothness-bound n)
        factor-base (cons -1 (primes-to (inc (int B))))
        relations   (add-relations n factor-base (inc (count factor-base)))]
    (find-factor n factor-base relations)))
