(ns numthy.arithmetic-fns
  (:require [clojure.math.combinatorics :refer [subsets]]
            [numthy.helpers :refer [isqrt divisible?]]
            [numthy.factorization.core :refer [distinct-prime-factors prime-factors-with-multiplicity]]
            [numthy.factorization.properties :refer [square-free?]]
            [numthy.primes.is-prime :refer [prime?]]))

(defn divisors-by-trial-division
  "Generates all proper divisors of a positive integer n (that is, all numbers < n
  that divide n without remainder) through trial division. Faster for small n (below ~500M)."
  [n]
  (if (= n 1) #{1}
    (->> (range 1 (inc (isqrt n)))
         (filter #(divisible? n %))
         (mapcat (fn [x] [x (/ n x)]))
         (remove #{n})
         (into (sorted-set)))))

(defn divisors-by-pollard
  "Generates all proper divisors of a positive integer n from its prime factorization.
  Faster for large n (above ~500M)."
  [n]
  (if (= n 1) #{1}
    (->> (prime-factors-with-multiplicity n)
         (subsets)
         (map #(reduce *' %))
         (drop-last)
         (into (sorted-set)))))

(defn divisors
  [n]
  (if (> n 500000000)
    (divisors-by-pollard n)
    (divisors-by-trial-division n)))

;; TODO: tau(n) returns the number of divisors of n
;; TODO: sigma(n) returns the sum of the divisors of n
;; TODO: sigmak(k,n) returns the sum of the kth power of the divisors of n

(defn chebyshev [x]
  ;; Second Chebyshev ψ(x) https://en.wikipedia.org/wiki/Chebyshev_function
  (->> (filter prime? (range (inc x)))
       (map #(*' (Math/floor (/ (Math/log x)
                                (Math/log %)))
                 (Math/log10 %)))
       (reduce +')))

(defn von-mangoldt
  [n]
  ;; Λ(n) https://en.wikipedia.org/wiki/Von_Mangoldt_function
  (let [pfs (distinct-prime-factors n)]
    (if (= 1 (count pfs))
      (Math/log10 (first pfs))
      0)))

(defn mobius
  [n]
  (if (square-free? n)
    (if (even? (count (distinct-prime-factors n))) 1 -1)
    0))
