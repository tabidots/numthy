(ns numthy.factorization.exponent-vectors
  (:require [clojure.math.numeric-tower :refer [expt]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.primes.sieves :as sv]
            [numthy.factorization.pollard :refer [brent-factorize]]))

(defn exp-vec-integer
  "Returns the exponent vector for the prime power representation of an integer,
  e.g., 168 = 2*2*2*3*7 = 2^3 * 3^1 * 5^0 * 7^1 -> [3 1 0 1]"
  [n]
  (when (and (integer? n) (pos? n))
    (if (= n 1) [0]
      (let [pfs        (brent-factorize n)
            largest-pf (-> pfs keys sort last)]
        (->> (sv/primes-to (inc largest-pf))
             (mapv #(get pfs % 0)))))))

(defn exp-vec-rational
  [q]
  (let [r (rationalize q)]
    ;; Sanity check to accept decimal representations of rational numbers
    ;; while still rejecting irrational numbers
    (when (or (ratio? q)
              (= (double r) (double q)))
      (let [n          (brent-factorize (numerator r))
            d          (->> (brent-factorize (denominator r))
                            (reduce-kv (fn [r k v] (assoc r k (- v))) {}))
            pfs        (merge-with + d n)
            largest-pf (-> pfs keys sort last)]
        (->> (sv/primes-to (inc largest-pf))
             (mapv #(get pfs % 0)))))))

(defn exponent-vector
  [n]
  (when (pos? n)
    (if (integer? n)
      (exp-vec-integer n)
      (exp-vec-rational n))))

(defn prime-powers->num
  "Converts an exponent vector into its corresponding integer or rational number."
  [pp]
  (let [primes (filter prime? (range))]
    (reduce *' (map expt primes pp))))
