(ns numthy.polynomials.core
  (:require [clojure.string :refer [join]]
            [numthy.modular-arithmetic.utils :refer [mod-inverse]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.arithmetic-fns :refer [divisors]]
            [clojure.math.numeric-tower :refer [abs lcm gcd]]))

(comment
  "To use these functions, represent polynomials as maps with powers as keys
  and coefficients as the corresponding vals. Terms with a coefficient of zero
  can be omitted (they will be trimmed anyway). Terms do not have to be in any
  particular order."
  "x^5 + 4x^3 + 27 ⟶" {0 27, 3 4, 5 1})

;; BASIC UTILITY FUNCTIONS

(defn- poly-trim
  "Removes terms with zero coefficients from a polynomial."
  [pnml]
  (->> (filter #(zero? (get pnml %)) (keys pnml))
       (reduce #(dissoc %1 %2) pnml)))

(defn degree
  "Finds the degree (power of highest-power term) of a polynomial."
  [pnml]
  (let [p (poly-trim pnml)]
    (if (empty? p) ##-Inf ;; the degree of the zero polynomial is −∞
      (apply max (keys p)))))

(defn lc
  "Leading coefficient (coefficient of highest-power term) of a polynomial."
  [pnml]
  (get pnml (degree (poly-trim pnml))))

(defn monic?
  "Tests if a polynomial is monic (that is, its leading coefficient is 1)."
  [pnml]
  (= 1 (lc pnml)))

;; SYMBOLIC REPRESENTATION

(defn- print-term
  "Prints a conventional symbolic representation of a single term in a polynomial."
  [power coeff]
  (let [x (cond
            (zero? power) ""
            (= 1 power)   "x"
            :else         (str "x^" power))
        c (cond
            (zero? power) coeff
            (= 1 coeff)   ""
            (= -1 coeff)  "-"
            :else         coeff)]
    (str c x)))

(defn poly-print
  "Prints a conventional symbolic representation of a polynomial."
  [pnml]
  (let [trimmed (poly-trim pnml)]
    (if (empty? trimmed) "0"
      (->> trimmed
           (sort)
           (reverse)
           (map (partial apply print-term))
           (join " + ")))))

;; BASIC ARITHMETIC OPERATIONS

(defn add
  "Adds two polynomials."
  [p1 p2]
  (poly-trim (merge-with + p1 p2)))

(defn mul
  "Multiplies two polynomials."
  [p1 p2]
  (->> (for [powers1 (keys p1) powers2 (keys p2)]
         {(+ powers1 powers2)
          (*' (get p1 powers1) (get p2 powers2))})
       (reduce #(merge-with + %1 %2) {})
       poly-trim))

(defn sub
  "Subtracts polynomial p2(x) from polynomial p1(x)."
  [p1 p2]
  (let [neg-p2 (mul p2 {0 -1})]
    (poly-trim (add p1 neg-p2))))

(comment
  "(x^2 + 2x + 1) * (x + 1) =  x^3 + 3x^2 + 3x + 1"
  (mul {0 1, 1 2, 2 1} {0 1, 1 1}) "=>" {0 1, 1 3, 2 3, 3 1})

;; LONG DIVISION

(defn- poly-quot
  "The quotient of a polynomial p1(x) divided by another p2(x).
  Returns nil if p2 is of a higher degree than p1."
  [p1 p2]
  (let [d1 (degree p1) d2 (degree p2)]
    (when (>= d1 d2)
      {(- d1 d2)                ;; power = difference in degree
       (/ (lc p1) (lc p2))})))  ;; coeff = quotient of lc's

(defn long-div
  "Polynomial long division of p1(x) / p2(x).
  Returns nil if p2 is of a higher degree than p1."
  ;; http://www.math.ucla.edu/~radko/circles/lib/data/Handout-358-436.pdf
  [p1 p2]
  (when-let [q (poly-quot p1 p2)]                ;; sanity check
    (loop [qs q r (->> q (mul p2) (sub p1))]
      (if (empty? r) {:quotient qs :remainder r} ;; divides evenly
        (if-let [new-q (poly-quot r p2)]
          (recur (conj qs new-q)                 ;; divides with remainder
                 (->> new-q (mul p2) (sub r)))
          {:quotient qs :remainder r})))))       ;; can't divide anymore

(defn div
  "Quotient of P(X) divided either by a rational divisor q or a polynomial Q(X)."
  [p q]
  (if (map? q) (long-div p q)
    (zipmap (keys p) (map #(/ % q) (vals p)))))

(comment
  "(x^2 - 4x + 4) / (x - 1) = (x - 3) with remainder 1."
  (div {0 4, 1 -4, 2 1} {0 -1, 1 1})
  "=>" {:quotient {1 1, 0 -3}, :remainder {0 1}})

(defn poly-rem
  "Remainder after dividing two polynomials, p1(x) / p2(x)."
  [p1 p2]
  (:remainder (div p1 p2)))

(defn quick-poly-rem
  "Shortcut to finding the remainder of p(x) / (x^r - 1).
  Slightly faster than doing it the normal way, but not by much."
  ;; Does this have something to do with "cyclotomic polynomials"?
  ;; https://en.wikipedia.org/wiki/Cyclotomic_polynomial
  ;; Source: https://medium.com/@sibu.it13/aks-primality-test-f184cf6365a1
  [pnml r]
  (poly-trim (reduce-kv (fn [res power coeff]
                           (merge-with + res {(mod power r) coeff}))
                        {} pnml)))

(comment
 "Compare times"
 (time (poly-rem (exp {0 1, 1 1} 5) {0 -1, 2 1}))
 (time (quick-poly-rem (exp {0 1, 1 1} 5) 2)))

;; MODULAR ARITHMETIC

(defn poly-mod
  "Reduces the terms of a given polynomial mod m."
  [pnml m]
  (zipmap (keys pnml) (map #(mod % m) (vals pnml))))

(defn exp
  "Exponentiation of a polynomial, [p(x)]^e."
  [pnml e]
  (reduce mul (repeat e pnml)))

(defn mod-exp
  "Slightly faster version of [p(x)]^e (mod m), where p(x) is a polynomial.
  Reduces the result of each multiplication mod m with every iteration, rather
  than only once at the end, in order to keep the intermediate coefficients
  from exploding."
  [pnml e m]
  (reduce #(poly-mod (mul %1 %2) m) (repeat e pnml)))

;; EVALUATION & INTERPOLATION

(defn poly-eval
  "Evaluates the value of a given polynomial at a point x using Horner's method."
  [pnml x]
  (let [powers (reverse (range (inc (degree pnml))))]
    (reduce (fn [sum p]
              (let [coeff (get pnml p 0)]
                (-> sum (* x) (+ coeff))))
            0 powers)))

(defn interpolate
  "Given a series of points as a hash-map, with x-coords as keys and y-coords
  as vals (e.g.,{0 -250, 10 0, 20 50, 30 -100}), uses Lagrange's method
  to interpolate the polynomial. Returns nil if there are duplicate x-coords.
  Will return integers or Ratios. Use mod-interpolate if you want only integers."
  [points]
  (let [xs (keys points)
        ys (vals points)]
    (when (distinct? xs)
      ;; http://wmueller.com/precalculus/families/lagrange.html
      ;; Transducer "basis-pnml" generate basis polynomials ℓ_0 ⋯ ℓ_k for each x-y pair.
      ;; Top is the product of polynomials (X - x_j) where X is the symbolic X (i.e., {1 1})
      ;; and x_j is each input x point (except x_i, the current).
      ;; Bottom is the product of (x_i - x_j) for all x_j except x_i (so none of the terms
      ;; are zero. In imperative terms, x_i is like the "outer loop" and x_j is like the
      ;; "inner loop".
      (letfn [(basis-pnml [xi]
                (let [top (reduce (fn [res xj]
                                    (if (= xi xj) res
                                      (mul res {1 1, 0 (- xj)})))
                                  {0 1} xs)
                      bot (reduce (fn [res xj]
                                    (if (= xi xj) res
                                      (*' res (- xi xj))))
                                  1 xs)]
                  (div top bot)))]
        (->> (map basis-pnml xs)       ;; Generate basis polynomials ℓ_j for each x,y pair
             (map #(mul {0 %1} %2) ys) ;; Scale each basis polynomial ℓ_j by constant y_j}
             (reduce add))))))         ;; Add the scaled basis polynomials

(defn mod-interpolate
  "Given a series of points as a hash-map, with x-coords as keys and y-coords
  as vals, and a modulus m (e.g.,{1 5, 2 6, 7 7} mod 2606193617), uses a
  modular-arithmetic version of Lagrange's method to interpolate a polynomial with
  only integer coefficients mod m. Returns nil if there are duplicate x-coords."
  [points m]
  (let [xs (keys points)
        ys (vals points)]
    (when (distinct? xs)
      (letfn [(mod-basis-pnml [xi]
                ;; https://www.judosaltgenius.com/2019/02/mod-squad
                (reduce (fn [res xj]
                          (if (= xi xj) res
                            (let [xixj-inv (mod-inverse (- xi xj) m)]
                              (mul res {1 xixj-inv,
                                        0 (-> xixj-inv (* (- xj)))}))))
                       {0 1} xs))]
        (let [raw (->> (map mod-basis-pnml xs)
                       (map #(mul {0 %1} %2) ys)
                       (reduce add))]
          (poly-mod raw m))))))

(defn mod-eval
  "Evaluates the value of a given polynomial at a point x, mod m, using Horner's method."
  [pnml x m]
  (-> (poly-eval pnml x) (mod m)))

;; CALCULUS

(defn differentiate
  "Given P(X), uses the power rule to find P'(X),"
  [pnml]
  (reduce-kv (fn [res power coeff]
               (if (zero? power) res
                 (merge-with + res {(dec power)
                                    (*' coeff power)})))
             {} (poly-trim pnml)))

;; ROOT-FINDING

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

;; TODO: Does divisors generate the same results as factors?
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

;; FACTORIZATION

(defn content
  "Returns the content of a polynomial (gcd of its coefficients) with either
  integer or rational coefficients."
  [pnml]
  (let [coefs (vals pnml)]
    (if (some ratio? coefs)
      ;; ↓ rational coefficients: LCM of denominators of all coefficients
      (let [d (reduce (fn [a b]
                        (lcm a (if (ratio? b) (denominator b) b)))
                      1 coefs)
            Q (mul pnml {0 d})]
        (/ (content Q) d))
      ;; ↓ all integer coefficients: GCD of coefficients
      (let [sign (if (neg? (lc pnml)) -1 1)]
        (* sign (reduce gcd coefs))))))

(defn ppc
  "Returns the primitive part–content factorization of a polynomial."
  [pnml]
  (let [c (content pnml) pp (div pnml c)]
    {:primitive-part pp :content c}))

;; TODO: https://en.wikipedia.org/wiki/Factorization_of_polynomials#Modern_methods
;; TODO: (irreducible? pnml)
