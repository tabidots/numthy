(ns numthy.factorization.dixon
  "Implementation of Dixon's integer factorization algorithm, which belongs to the family of congruence-of-squares methods."
  (:require [clojure.core.reducers :as r]
            [numthy.helpers :refer [isqrt]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.factorization.squares-utils :refer [smooth? smoothness-bound find-factor-from-congruent-relations]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

(defn- generate-relations
  "Given a factor base (vector of primes) with B elements, accumulate B+1 values of z,
  between sqrt(n) and n, s.t. z^2 mod n is smooth over the factor base. Returns a hash-map
  with the zs as keys and the exponent vectors as vals."
  [n factor-base]
  (let [prime-map   (into (sorted-map) (zipmap factor-base (repeat 0)))
        z2-Bsmooth? (fn [z]
                      (when-let [expv (some-> (mod-pow z 2 n)
                                              (smooth? factor-base prime-map))]
                        {z expv}))
        relations   (->> (range (inc (isqrt n)) (inc n))
                         (keep z2-Bsmooth?)
                         (take (inc (count factor-base))))]
    (apply merge relations)))

(defn dixon-factorize
  "Uses Dixon's factorization algorithm to find at least one factor of n.
  Returns nil if algorithm fails to find a factor."
  ;; With reference to https://crypto.stanford.edu/cs359c/17sp/projects/BrendonGo.pdf
  [n]
  (let [B           (smoothness-bound n)
        factor-base (filter prime? (range B))
        relations   (generate-relations n factor-base)]
    (find-factor-from-congruent-relations n factor-base relations)))
