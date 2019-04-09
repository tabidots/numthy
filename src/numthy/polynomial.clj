(ns numthy.polynomial
  (:require [clojure.string :refer [join]]
            [numthy.modular-arithmetic :refer [mod-inverse]]))

(comment
  "To use these functions, represent polynomials as maps with powers as keys
  and coefficients as the corresponding vals. Terms with a coefficient of zero
  can be omitted (they will be trimmed anyway). Terms do not have to be in any
  particular order."
  "x^5 + 4x^3 + 27 ⟶" {0 27, 3 4, 5 1})

;; BASIC UTILITY FUNCTIONS

(defn poly-trim-
  "Removes terms with zero coefficients from a polynomial."
  [pnml]
  (->> (filter #(zero? (get pnml %)) (keys pnml))
       (reduce #(dissoc %1 %2) pnml)))

(defn degree
  "Finds the degree (power of highest-power term) of a polynomial."
  [pnml]
  (apply max (keys pnml)))

(defn lc
  "Leading coefficient (coefficient of highest-power term) of a polynomial."
  [pnml]
  (get pnml (degree pnml)))

(defn monic?
  "Tests if a polynomial is monic (that is, its leading coefficient is 1)."
  [pnml]
  (= 1 (lc pnml)))

;; SYMBOLIC REPRESENTATION

(defn print-term-
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
  (let [trimmed (poly-trim- pnml)]
    (if (empty? trimmed) "0"
      (->> trimmed
           (sort)
           (reverse)
           (map (partial apply print-term-))
           (join " + ")))))

;; BASIC ARITHMETIC OPERATIONS

(defn add
  "Adds two polynomials."
  [p1 p2]
  (poly-trim- (merge-with + p1 p2)))

(defn mul
  "Multiplies two polynomials."
  [p1 p2]
  (->> (for [powers1 (keys p1) powers2 (keys p2)]
         {(+ powers1 powers2)
          (*' (get p1 powers1) (get p2 powers2))})
       (reduce #(merge-with + %1 %2) {})
       poly-trim-))

(defn sub
  "Subtracts polynomial p2(x) from polynomial p1(x)."
  [p1 p2]
  (let [neg-p2 (mul p2 {0 -1})]
    (poly-trim- (add p1 neg-p2))))

(comment
  "(x^2 + 2x + 1) * (x + 1) =  x^3 + 3x^2 + 3x + 1"
  (mul {0 1, 1 2, 2 1} {0 1, 1 1}) "=>" {0 1, 1 3, 2 3, 3 1})

;; LONG DIVISION

(defn poly-quot-
  "The quotient of a polynomial p1(x) divided by another p2(x).
  Returns nil if p2 is of a higher degree than p1."
  [p1 p2]
  (let [d1 (degree p1) d2 (degree p2)]
    (when (>= d1 d2)
      {(- d1 d2)                ;; power = difference in degree
       (/ (lc p1) (lc p2))})))  ;; coeff = quotient of lc's

(defn div
  "Polynomial long division of p1(x) / p2(x).
  Returns nil if p2 is of a higher degree than p1."
  ;; http://www.math.ucla.edu/~radko/circles/lib/data/Handout-358-436.pdf
  [p1 p2]
  (when-let [q (poly-quot- p1 p2)]               ;; sanity check
    (loop [qs q r (->> q (mul p2) (sub p1))]
      (if (empty? r) {:quotient qs :remainder r} ;; divides evenly
        (if-let [new-q (poly-quot- r p2)]
          (recur (conj qs new-q)                 ;; divides with remainder
                 (->> new-q (mul p2) (sub r)))
          {:quotient qs :remainder r})))))       ;; can't divide anymore

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
  (poly-trim- (reduce-kv (fn [res power coeff]
                           (merge-with + res {(mod power r) coeff}))
                         {} pnml)))

(comment
 (time (poly-rem (exp {0 1, 1 1} 5) {0 -1, 2 1}))
 "Elapsed time: 0.261586 msecs"
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
                  (:quotient (div top {0 bot}))))]
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
