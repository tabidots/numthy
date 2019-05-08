(ns numthy.factorization.trial-division
  "Batch trial division, gcd, and smooth-filtering functions based on algorithms by Dan Bernstein."
  (:require [clojure.math.numeric-tower :refer [abs expt gcd]]
            [clojure.core.reducers :as r]
            [numthy.modular-arithmetic.utils :refer [mod-pow]]))

(defn- product-tree
  "Iteratively consumes a collection of integers in pairs to find the product
  of the collection. For very large collections, this is much faster than (reduce * coll)."
  ;; https://facthacks.cr.yp.to/product.html
  [coll]
  (letfn [(twin-prods [c]
            (when (> (count c) 1)
              (mapv (fn [[a b]] (*' a b)) (partition 2 2 [1] c))))]
    (vec (take-while some? (iterate twin-prods coll)))))

(defn- fast-product
  [coll]
  (first (peek (product-tree coll))))

(defn- remainder-tree
  "Walks down through the product tree of coll to find n mod x for every x in coll.
  For very large collections, this is faster than (map #(mod n %) coll)."
  ;; https://facthacks.cr.yp.to/remainder.html
  [n coll]
  (letfn [(twin-rems [as bs]
            (mapv #(mod %1 %2) (mapcat #(repeat 2 %) as) bs))]
    (->> (product-tree coll)
         reverse
         vec
         (r/reduce (fn [a b] (twin-rems a b)) [n]))))

(defn- naive-batch-gcd
  "Finds the GCD of every integer in coll with the product of all other integers in coll.
  Works well for small collections."
  [coll]
  (let [product    (r/fold *' coll)
        remainders (pmap #(mod product %) (pmap #(expt % 2) coll))]
    (pmap (fn [r x] (gcd (/ r x) x)) remainders coll)))

(defn batch-gcd
  "Finds the GCD of every integer in coll with the product of all other integers in coll.
  Massive speedup for very large collections."
  ;; https://facthacks.cr.yp.to/batchgcd.html
  [coll]
  (let [remainders (remainder-tree (fast-product coll) (mapv #(expt % 2) coll))]
     (pmap (fn [r x] (gcd (/ r x) x)) remainders coll)))

(defn cocomposites-in
  "Filters the integers in coll that have a shared factor with at least one other integer in coll."
  [coll]
  (->> (pmap (fn [x gcd'] (when (> gcd' 1) x)) coll (batch-gcd coll))
       (remove nil?)
       (into [])))

(defn- primes-in-product
  "Returns the subset of primes that divide the product of coll, in a vector. Note that
  there may be other primes that divide the product of coll."
  [primes coll]
  (let [remainders (remainder-tree (fast-product coll) primes)]
    (->> (map (fn [pr rm] (when (zero? rm) pr)) primes remainders)
         (remove nil?)
         (into []))))

(defn- primes-in-each
  "Recursively finds the subsets [s0, s1...] of primes [p0, p1...] that divide [x0, x1...]."
  ;; https://facthacks.cr.yp.to/batchtrial.html
  [primes coll]
  (if (empty? coll) []
    (let [primes' (primes-in-product primes coll)
          c       (count coll)]
      (if (= c 1) [primes']
        (mapcat #(primes-in-each primes' %) (split-at (quot c 2) coll))))))

(defn batch-trial-division
  "Recursively finds the subsets [s0, s1...] of primes [p0, p1...] that divide
  [x0, x1...] and returns them in a map of the form {x0 s0, x1 s1...}."
  ;; https://facthacks.cr.yp.to/batchtrial.html, http://cr.yp.to/papers/sf.pdf
  [primes coll]
  (zipmap coll (primes-in-each primes coll)))

(defn batch-smooth-filter
  "Given a large collection of candidates and a collection of primes, returns those
  candidates which factor completely over the collection of primes. Optionally, if
  a transducer is given, returns a hash-map where the values are the transduced candidates
  that factor completely and the keys are the original candidates that produced those values.
  This function allows transducers that return nil."
  ;; Adapted from https://hal.inria.fr/inria-00188645v1/document
  ([candidates primes]
   (let [remainders (remainder-tree (fast-product primes) candidates)]
     (remove nil?
       (pmap (fn [rm cand]
               (let [e (first (filter #(>= (expt 2 (* 2 %)) cand) (rest (range))))
                     ;; ↑ for each remainder z_i and candidate x_i, {z_i}^2^e where
                     ;; e is the smallest pos-int s.t. 2^2^e ≧ {x_i}
                     y (mod-pow rm (* 2 e) cand)]
                   (when (= cand (gcd cand y)))))
             remainders candidates))))
  ([raw-candidates xf primes]
   (let [[candidates x-candidates] ((juxt filter keep) xf raw-candidates)
         remainders (remainder-tree (fast-product primes) candidates)]
     (apply merge
       (pmap (fn [rm cand x-cand]
               (let [e (first (filter #(>= (expt 2 (* 2 %)) x-cand) (rest (range))))
                     y (mod-pow rm (* 2 e) x-cand)]
                (when (= x-cand (gcd x-cand y))
                  {cand x-cand})))
            remainders candidates x-candidates)))))
