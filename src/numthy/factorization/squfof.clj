(ns numthy.factorization.squfof
  "Implementation of Shanks's square forms factorization algorithm."
  (:require [clojure.math.numeric-tower :refer [expt gcd sqrt]]
            [numthy.helpers :refer [isqrt]]
            [numthy.perfect-powers :refer [perfect-square?]]
            [numthy.primes.sieves :refer [squarefrees-to]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

(declare squfof)

(defn squfof-single
  [n k]
  (let [kn                (*' k n)
        bound             (-> (expt n 1/5) (quot 2))
                          ;; ↑ Taken from Tilman Neumann, minus the `cmult` parameter
                          ;; Divide by 2 b/c only taking alternate iterations
        p0                (isqrt kn)
        q0                (- kn (*' p0 p0))]
     (if (zero? q0) (squfof p0)
       (let [step             (fn [[p-old p b q-cur q-old]]
                                (let [b      (-> (+' p0 p) (/ q-cur) bigint)
                                      p'     (-> (*' b q-cur) (- p))
                                      q-cur' (-> (- p p') (*' b) (+' q-old))]
                                  [p p' b q-cur' q-cur]))
             [p-old _ _ _ qi] (->> [nil p0 nil q0 1]
                                   (iterate step)
                                   (take-nth 2)      ;; take every even iteration
                                   rest              ;; drop the 0th iteration
                                   (take bound)
                                   (filter #(perfect-square? (peek %)))
                                   first)]           ;; and go until q-cur is a square
         (when (some? qi)
           (let [sqrt-qi (isqrt qi)
                 b0      (-> (- p0 p-old) (/ sqrt-qi) bigint)
                 p0      (-> (*' b0 sqrt-qi) (+' p-old))
                 q1      (-> (- kn (*' p0 p0)) (/ sqrt-qi))
                 p       (->> (iterate step [nil p0 b0 q1 sqrt-qi])
                              (filter #(= (first %) (second %)))
                              ffirst)
                 factor  (gcd n p)]
             (when (< 1 factor n) factor)))))))

(def multipliers
  "Based on Gower & Wagstaff's finding that the optimal multiplier k for any n
  to be factored is among the divisors of 1155, the product of the first four odd primes."
  ;; https://homes.cerias.purdue.edu/~ssw/squfof.pdf
  [1 3 5 7 11 15 21 33 35 55 77 105 165 231 385 1155])

(def second-chance
  ;; Sometimes, the divisors of 1155 don't find a factor! e.g. 11529083
  ;; Even 100 fails sometimes (squfof finds a factor of 6518056469576923 when k = 286, 533...)
  (remove (set multipliers) (squarefrees-to 100)))

(defn squfof
  [n]
  (when (>= n 1000)
    (or (some #(squfof-single n %) multipliers)
        (some #(squfof-single n %) second-chance))))
