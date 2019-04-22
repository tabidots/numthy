(ns numthy.polynomials.cyclotomic
  (:require [numthy.polynomials.core :refer [mul div]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.arithmetic-fns :refer [divisors mobius]]))

;; CYCLOTOMIC POLYNOMIALS
;; TODO: Figure out what these actually are and what can be done with them
;; TODO: Does factors output the same result as divisors?

(defn cyclotomic
  "Returns the nth cyclotomic polynomial, i.e., the unique irreducible polynomial
  with integer coefficients that is a divisor of x^n - 1 and is not a divisor of
  x^k - 1 for any k < n. This function uses the Möbius inversion formula to
  generate the polynomial."
  [n]
  (when (pos? n)
    (if (= 1 n) {0 -1, 1 1}
      (let [{top-roots 1 bot-roots -1} (group-by mobius (divisors n))
            top (reduce mul (map (fn [d] {(/ n d) -1, 0 1}) top-roots))
            bot (reduce mul (map (fn [d] {(/ n d) -1, 0 1}) bot-roots))]
        (:quotient (div top bot))))))

(defn palindromic?
  "A polynomial is palindromic if its coefficients, in order of powers of the terms,
  are the same forward and backward."
  [pnml]
  (let [coefs (vals (into (sorted-map) pnml))]
    (= coefs (reverse coefs))))

(defn anti-palindromic?
  "A polynomial is anti-palindromic if its coefficients, in order of powers of
  the terms, are the negatives of themselves when the order is reversed."
  [pnml]
  (let [coefs (vals (into (sorted-map) pnml))]
    (= coefs (map - (reverse coefs)))))
