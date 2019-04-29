(ns numthy.factorization.brent
  (:require [clojure.math.numeric-tower :refer [gcd expt abs]]
            [numthy.helpers :refer [rand-num]]
            [numthy.primes.is-prime :refer [quick-prime?]]
            [numthy.modular-arithmetic.utils :refer [mod-pow mod-mul]]
            [clojure.core.reducers :as r]))

(defn- pollard-brent
  "Uses a slightly abbreviated version of the Pollard-Brent algorithm to find one
  prime factor of a large composite integer n. This leverages the optimization of
  r/fold to dispense with the need to have a step size `m` and to keep and calculate
  all intermediate values of `q` based on `m`."
  [n]
  (if (even? n) 2
    (let [const (rand-num 1 (dec n))
          step  #(-> (+ const (mod-pow % 2 n)) (mod n))
          seed  (rand-num 1 (dec n))]
      (loop [distance 1 tortoise seed q 1]
        (let [hare      (nth (iterate step tortoise) distance)
              hare-path (into [] (r/take distance (iterate step (step hare))))
              q*        (->> hare-path
                             (r/map #(abs ^BigInteger (- tortoise %)))
                             (r/fold (fn ([] q)
                                         ([q z] (mod-mul q z n)))))
              factor   (gcd ^BigInteger q* n)]
          (cond
            (< 1 factor n) factor
            (= factor n)   (->> hare-path  ;; Backtracking
                                (r/map #(gcd ^BigInteger (abs ^BigInteger (-' tortoise %)) n))
                                (r/filter #(< 1 % n))
                                (r/foldcat)
                                first)
            ;; If distance gets too high, we are dealing with very large prime factors.
            ;; (> distance 4000000)  (some-other-factorization-method n)
            :else          (recur (*' distance 2) (peek hare-path) q*)))))))

(defn brent-factorize
  "Uses the Pollard-Brent algorithm repeatedly to find all prime factors of a large
  integer n, along with their powers, in a hash-map of the form {p_0 e_0, ⋯, p_i e_i}."
  [n]
  (when (and (integer? n) (> n 1))
    (loop [n n res (sorted-map)]
      (if (quick-prime? n)
        (merge-with + res {n 1})
        (let [p (first (filter quick-prime? (repeatedly #(pollard-brent n))))]
          (recur (/ n p) (merge-with + res {p 1})))))))
