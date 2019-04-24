(ns numthy.factorization.core
  (:require [clojure.math.numeric-tower :refer [gcd expt]]
            [numthy.helpers :refer [rand-num]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]))

;; TODO: Quadratic sieve algorithm

(defn- pollard-rho
  "Uses Pollard's rho algorithm to find one prime factor of a large composite integer n.
  Terminates only if n is composite."
  [n]
  (if (even? n) 2 ;; Pollard's rho can only return odd numbers
    (let [r (rand-num 1 (dec n))
          g (fn [x] (+ r (mod-pow x 2 n)))]    ;; g(x) = (x^2 + random) mod n
      (loop [a (g 2) b (g (g 2)) d 1]
        (cond
          (= d n)   (pollard-rho n) ;; keep trying random numbers till a factor is found
          (< 1 d n) d
          :else     (recur
                      (g a) (g (g b)) (gcd (- a b) n)))))))

(defn pollard-factorize
  "Uses Pollard's rho algorithm recursively to find all prime factors of a large
  integer n, along with their powers, in a hash-map of the form {p_0 e_0, ⋯, p_i e_i}."
  [n]
  (when (and (integer? n) (> n 1))
    (loop [n n res (sorted-map)]
      (if (.isProbablePrime (biginteger n) 5)
        (merge-with + res {n 1})
        (let [p (first (filter prime? (repeatedly #(pollard-rho n))))]
          (recur (/ n p) (merge-with + res {p 1})))))))

(defn distinct-prime-factors
  "Returns the distinct prime factors of an integer n. e.g., 1300 = (2 5 13)"
  [n]
  (when-let [pfact (pollard-factorize n)]
    (keys pfact)))

(defn prime-factors-with-multiplicity
  "Returns a vector of the prime factors of an integer n with multiplicity.
  e.g., 168 = [2 2 2 2 3 7]"
  [n]
  (reduce-kv (fn [r b e]
               (into r (repeat e b)))
             [] (pollard-factorize n)))

(defn prime-signature
  "Returns the non-zero exponents of the prime factorization of n."
  [n]
  (when-let [pfact (pollard-factorize n)]
    (vals pfact)))

;; TODO: Implement Brent's algorithm
;; https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf

(defn phi
  "Euler's totient function, optimized using the product rule and Pollard's rho algorithm
  for super fast factorization."
  [n]
  (->> (distinct-prime-factors n)
       (map #(- 1 (/ 1 %)))
       (reduce * n)))
