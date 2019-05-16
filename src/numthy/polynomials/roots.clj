(ns numthy.polynomials.roots
  (:require [numthy.arithmetic-fns :refer [divisors]]
            [numthy.polynomials.core :refer [degree lc poly-rem div]]
            [clojure.math.numeric-tower :refer [abs]]))

(defn- quadratic-roots
  "Given a quadratic polynomial, use the quadratic formula to find its
  real-valued zeroes (nil if none)."
  [pnml]
  (when (= 2 (degree pnml))
    (let [a (get pnml 2 0) b (get pnml 1 0) c (get pnml 0 0)
          discriminant (-> (*' b b)
                           (- (* 4 a c)))]
      (cond
        (pos? discriminant)  [(-> (- b) (+ (Math/sqrt discriminant))
                                  (/ (* 2 a)))
                              (-> (- b) (- (Math/sqrt discriminant))
                                  (/ (* 2 a)))]
        (zero? discriminant) [(-> (- b) (/ (* 2 a)))]
        :else nil))))

(defn- possible-roots
  "Uses the Rational Root Theorem to generate all possible rational roots of a polynomial."
  [pnml]
  (if-let [constant (get pnml 0)]
    (->> (for [p (divisors (abs constant))
               q (divisors (abs (lc pnml)))]
           [(/ p q) (- (/ p q))])
         (apply concat)
         (into (sorted-set)))
    #{0}))

(defn- root?
  "r is a root of a polynomial if no remainder is left after dividing pnml / (x - r)."
  [pnml r]
  (empty? (poly-rem pnml {1 1, 0 (- r)})))

(defn roots
  "Finds all rational and quadratic roots of a polynomial of arbitrary degree.
  Cannot return more than two irrational roots. Returns nil if there are no
  such roots, or if the polynomial is of degree 0."
  ;; https://courses.lumenlearning.com/waymakercollegealgebra/chapter/zeros-of-a-polynomial-function/
  [pnml]
  (loop [pnml  pnml
         roots (sorted-set)]
    (condp = (degree pnml)
      0 nil
      1 (conj roots (/ (get pnml 0 0) (- (get pnml 1))))
      2 (into roots (quadratic-roots pnml))
      (let [winners (filter #(root? pnml %) (possible-roots pnml))]
        (if (empty? winners) ;; No solutions for the current iterations
          (if (empty? roots) nil roots) ;; Check whether solutions at all
                 ;; ↓ sucessively divide pnml by (x - r)
          (recur (:quotient (div pnml {1 1, 0 (- (first winners))}))
                 (into roots winners)))))))
