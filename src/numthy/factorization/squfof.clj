(ns numthy.factorization.squfof
  "Implementation of Shanks's square forms factorization algorithm."
  (:require [clojure.math.numeric-tower :refer [expt gcd sqrt]]
            [numthy.helpers :refer [isqrt]]
            [numthy.perfect-powers :refer [perfect-square?]]
            [numthy.primes.sieves :refer [squarefrees-to]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

(defn squfof-single
  [n k]
  (let [kn                (*' k n)
        bound             (-> (expt n 1/5) (quot 2))
                          ;; ↑ Taken from Tilman Neumann, minus the `cmult` parameter
                          ;; Divide by 2 b/c only taking alternate iterations
        p0                (isqrt kn)
        step              (fn [[p-old p b q-cur q-old]]
                            (let [b      (-> (+' p0 p) (/ q-cur) bigint)
                                  p'     (-> (*' b q-cur) (- p))
                                  q-cur' (-> (- p p') (*' b) (+' q-old))]
                              [p p' b q-cur' q-cur]))
        [p-old _ _ _ qi]  (->> [nil p0 nil (- kn (*' p0 p0)) 1]
                               (iterate step)
                               (take-nth 2)      ;; take every even iteration
                               rest              ;; drop the 0th iteration
                               (take bound)
                               (filter #(perfect-square? (peek %)))
                               first)]           ;; and go until q-cur is a square
    (when (some? qi)
      (let [sqrt-qi           (isqrt qi)
            b0                (-> (- p0 p-old) (/ sqrt-qi) bigint)
            p0                (-> (*' b0 sqrt-qi) (+' p-old))
            q1                (-> (- kn (*' p0 p0)) (/ sqrt-qi))
            p                 (->> (iterate step [nil p0 b0 q1 sqrt-qi])
                                   (filter #(= (first %) (second %)))
                                   ffirst)
            factor            (gcd n p)]
        (when (< 1 factor n) factor)))))

(defn squfof
  [n]
  ;; Don't bother with generating square-free multipliers b/c it takes time to calculate them
  ;; Not really clear how many multipliers to try
  (some #(squfof-single n %) (rest (range 100))))
