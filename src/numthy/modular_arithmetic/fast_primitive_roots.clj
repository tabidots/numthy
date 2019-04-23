(ns numthy.modular-arithmetic.fast-primitive-roots
  (:require [clojure.math.numeric-tower :refer [expt]]
            [clojure.core.reducers :as r]
            [numthy.helpers :refer [coprime? rand-num]]
            [numthy.modular-arithmetic.utils :refer [mod-pow odd-prime? odd-prime-power?]]
            [numthy.modular-arithmetic.groups :refer [multiplicative-group]]
            [numthy.factorization.core :as f]))

;; Rewriting the original primitive roots functions using clojure.core.reducers
;; Significant speedup, but only if the result is left as an unordered vector

(defn lpr-prime
  [p]
  (if (= p 2) 1
    (when (odd-prime? p)
      (let [phi     (f/phi p);(dec p)
            pfs     (f/distinct-prime-factors phi)
            p-root? (fn [a]
                      (and (not-any? #(= 1 (mod-pow a (/ phi %) p)) pfs)))]
        (first (filter p-root? (drop 2 (range))))))))

(defn r-proots-prime
  [p]
  (let [lpr (lpr-prime p) phi (f/phi p)]
    (->> (range 1 phi)
         (r/filter #(coprime? % phi))
         (r/map #(mod-pow lpr % p))
         (r/foldcat))))

(defn- r-proots-prime-squared
  [p]
  (let [n    (*' p p)
        lift (fn [a] (r/map #(+ a (* % p)) (range p)))]
    (->> (r-proots-prime p)
         (r/mapcat (fn [r]
                     (->> (lift r)
                          (r/remove #(= 1 (mod-pow % (dec p) n))))))
         (r/foldcat))))

(defn r-proots-prime-power
  [[p k]]
  (cond
    (= k 1)                 (r-proots-prime p)
    (and (= p 2) (> k 2))   nil
    :else
    (letfn [(lift [a e] (r/map #(+ a (* % (expt p e))) (range p)))]
      (loop [j 2 roots (r-proots-prime-squared p)]
        (if (= j k) roots
          (recur (inc j)
                 (->> roots
                      (r/mapcat #(lift % j))
                      (r/foldcat))))))))

(defn r-proots-2pk
  [[p k]]
  (when-let [roots (r-proots-prime-power [p k])]
    (let [pk (expt p k)
          n  (* 2 pk)]
      (->> roots
           (r/map #(if (odd? %) % (mod (+ % pk) n)))
           (r/foldcat)))))

(defn r-primitive-roots
  [n]
  (cond
    (< n 1) nil
    (= n 1) #{0}
    (= n 2) #{1}
    (= n 4) #{3}
    :else
    (let [opp       (odd-prime-power? n)
          twice-opp (odd-prime-power? (* n 1/2))]
      (cond
        (some? opp)       (r-proots-prime-power opp)
        (some? twice-opp) (r-proots-2pk twice-opp)
        :else             nil))))
