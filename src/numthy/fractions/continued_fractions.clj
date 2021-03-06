(ns numthy.fractions.continued-fractions
  (:require [clojure.math.numeric-tower :refer [sqrt]]
            [numthy.helpers :refer [isqrt]]))

(defn rational-cfrac
  "Continued fraction representation of n, where n is a rational number.
  This can work for irrational numbers too, but seems to lose accuracy
  after a certain number of iterations."
  ; https://en.wikipedia.org/wiki/Continued_fraction
  [n]
  (when (rational? n)
    (letfn [(cfrac [m]
              (loop [a (rationalize m) cf []]
                (if (not (ratio? a)) ;; finite continued fraction
                  (conj cf a)
                  (let [q (quot (numerator a) (denominator a))]
                    (recur (/ 1 (rem a q))
                           (conj cf q))))))]
      (if (< 0 n 1)
        (cons 0 (map int (cfrac (/ 1 n))))
        (map int (cfrac n))))))

(defn sqrt-cfrac
  "Continued fraction representation of sqrt(n) where sqrt(n) is irrational.
  Outputs the starting term and one complete period. Works well up to int limit."
  ; http://www.maths.surrey.ac.uk/hosted-sites/R.Knott/Fibonacci/cfINTRO.html#section6.2.2
  [n]
  (when (and (int? n)
             (not (rational? (sqrt n))))           ; sanity check
    (let [a (isqrt n)]
      (loop [top 0 bot 1 cf []]
        (let [q         (quot (+ a top) bot)
              conjugate (- (- top (* q bot)))
              denom     (/ (- n (* conjugate conjugate))
                           bot)]   ; divide by bot to reduce the denom to its lowest term
          (if (= denom 1)
            (conj cf q (+ a conjugate))
            (recur conjugate denom (conj cf q))))))))

(defn lazy-sqrt-cfrac
  [n]
  (let [a0 (isqrt n)
        cfrac-next (fn [x]
                     (when-some [[m d a] x]
                       (when (or (zero? m) (> d 1))
                         (let [m' (-> (* d a) (- m))
                               d' (/ (- n (* m' m')) d)
                               a' (quot (+ a0 m') d')]
                           [m' d' a']))))]
    (lazy-seq (map peek (iterate cfrac-next [0 1 a0])))))

(defn ->number
  "Converts a continued fraction representation to an integer or Ratio."
  [cf]
  (->> (reverse cf)
       (reduce (fn [prev-sum t] (+ t (/ 1 prev-sum))))))

(defn c-step
  [x [a b]]
  (+ (* x a) b))

(defn ->convergents
  "Returns the rational convergents of a given continued fraction."
  [cfc]
  (loop [as cfc hs '(1 0) ks '(0 1)] ;; Reverse order because we are cons'ing
    (if-let [a (first as)]
      (recur (rest as)
        (cons (c-step a (take 2 hs)) hs)
        (cons (c-step a (take 2 ks)) ks))
      (map / (drop 2 (reverse hs))
           (drop 2 (reverse ks))))))

;; TODO: More stuff from http://www.maths.surrey.ac.uk/hosted-sites/R.Knott/Fibonacci/cfINTRO.html
;; http://www.maths.surrey.ac.uk/hosted-sites/R.Knott/Fibonacci/cfCALCbn.html

;; https://en.wikipedia.org/wiki/Periodic_continued_fraction
;; https://en.wikipedia.org/wiki/Quadratic_irrational_number
;; https://en.wikipedia.org/wiki/Solving_quadratic_equations_with_continued_fractions
;; https://en.wikipedia.org/wiki/Generalized_continued_fraction#Roots_of_positive_numbers

(comment

 ;; https://en.wikipedia.org/wiki/Metallic_mean
 "In addition to the ratio (1 + sqrt(5)) / 2,
 the golden mean can be calculated to a high precision from its
 continued fraction representation, [2; 2, 2, 2, 2...]"

 (with-precision 1000 (* 1M (expand-cfrac (repeat 300 1)))))
