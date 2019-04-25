(ns numthy.fractions.continued-fractions
  (:require [clojure.math.numeric-tower :refer [sqrt]]))

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
  Outputs the starting term and one complete period."
  ; http://www.maths.surrey.ac.uk/hosted-sites/R.Knott/Fibonacci/cfINTRO.html#section6.2.2
  [n]
  (when (and (int? n)
             (not (rational? (sqrt n))))           ; sanity check
    (let [a (last (take-while #(<= (* % %) n) (range)))] ; a for "approximation"
      (loop [top 0 bot 1 cf []]
        (let [q         (quot (+ a top) bot)
              conjugate (- (- top (* q bot)))
              denom     (/ (- n (* conjugate conjugate))
                           bot)]   ; divide by bot to reduce the denom to its lowest term
          (if (= denom 1)
            (conj cf q (+ a conjugate))
            (recur conjugate denom (conj cf q))))))))

(defn cfrac->decimal
  "Converts a continued fraction representation to decimal."
  [cf]
  (->> (reverse cf)
       (reduce (fn [prev-sum t] (+ t (/ 1 prev-sum))))))

(comment

 ;; https://en.wikipedia.org/wiki/Metallic_mean
 "In addition to the ratio (1 + sqrt(5)) / 2,
 the golden mean can be calculated to a high precision from its
 continued fraction representation, [2; 2, 2, 2, 2...]"

 (with-precision 1000 (* 1M (expand-cfrac (repeat 300 1)))))
