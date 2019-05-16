(ns numthy.factorization.pollard
  (:require [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :refer [gcd expt abs]]
            [numthy.helpers :refer [rand-num]]
            [numthy.modular-arithmetic.utils :refer [mod-pow mod-mul]]
            [numthy.primes.is-prime :refer [quick-prime?]]
            [numthy.primes.sieves :refer [primes-to]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

(defn rho
  "Uses Pollard's rho algorithm to find one prime factor of a large composite integer n.
  Terminates only if n is composite. Returns nil if the algorithm fails."
  [n]
  (if (even? n) 2 ;; Pollard's rho can only return odd numbers
    (let [c (rand-num 1 (dec n))
          g (fn [x] (-> (+' c (mod-pow x 2 n)) (mod n)))]    ;; g(x) = (x^2 + random) mod n
      (loop [a (g 2) b (g (g 2)) d 1N]
        (cond
          (= d n)   nil ;; Failure
          (< 1 d n) d
          :else     (recur (g a) (g (g b)) (gcd n (-> (abs (-' a b))
                                                      (mod n)))))))))
(defn brent
  "Uses a slightly abbreviated version of the Pollard-Brent algorithm to find one
  prime factor of a large composite integer n. This leverages the optimization of
  r/fold to dispense with the need to have a step size `m` and to keep and calculate
  all intermediate values of `q` based on `m`."
  ([n]
   (brent n 4000000))
  ([n limit]
   (if (even? n) 2
     (let [const (rand-num 1 (max 2 (dec n)))
           step  #(-> (+ const (mod-pow % 2 n)) (mod n))
           seed  (rand-num 1 (max 2 (dec n)))]
       (loop [distance 1N tortoise seed q 1N]
         (let [hare            (nth (iterate step tortoise) distance)
               ;; Interesting: This implementation is not an exact recreation of Brent's algo but seems to work faster than the proper version
               hare-path       (into [] (r/take distance (drop 1 (iterate step hare))))
               ^BigInteger q*  (->> hare-path
                                    (r/map #(abs ^BigInteger (- tortoise %)))
                                    (r/fold (fn ([] q)
                                               ([q z] (mod-mul q z n)))))
               factor          (gcd ^BigInteger q* n)]
           (cond
             (Thread/interrupted) (throw (InterruptedException. "Interrupting Pollard-Brent..."))
             (> distance limit)   nil
             (< 1 factor n)       factor
             (= factor n)         (->> hare-path  ;; Backtracking
                                       (r/map #(gcd ^BigInteger (abs ^BigInteger (-' tortoise %)) n))
                                       (r/filter #(< 1 % n))
                                       (r/foldcat)
                                       first)
             :else                (recur (*' distance 2) (peek hare-path) q*))))))))

(defn fast-brent
  [n]
  (brent n 17000))

(defn brent-factorize
  "Uses the Pollard-Brent algorithm repeatedly to find all prime factors of a large
  integer n, along with their powers, in a hash-map of the form {p_0 e_0, ⋯, p_i e_i}."
  [n]
  (when (and (integer? n) (> n 1))
    (loop [n n res (sorted-map)]
      (if (quick-prime? n)
        (merge-with + res {n 1})
        (let [p (first (filter quick-prime? (repeatedly #(brent n))))]
          (recur (/ n p) (merge-with + res {p 1})))))))

(comment
 "Accelerated version, as used above. Not sure why it works faster, but it does."
 (loop [distance 1 tortoise 1]
   (let [_         (println "Distance:" distance)
         _         (println "Tortoise at:" tortoise)
         hare      (nth (iterate inc tortoise) distance)
         _         (println "Hare starts at:" hare)
         hare-path (into [] (r/take distance (drop 1 (iterate inc hare))))
         _         (println "Hare path:" hare-path)]
     (cond
       (> distance 16) nil
       :else           (recur (*' distance 2) (peek hare-path)))))

 "Correct version. Makes more sense at a glance, but is slower"
 (loop [distance 1 tortoise 1]
   (let [_         (println "Distance:" distance)
         _         (println "Tortoise at:" tortoise)
         hare-path (into [] (r/take distance (drop 1 (iterate inc tortoise))))
         _         (println "Hare path:" hare-path)]
     (cond
       (> distance 16) nil
       :else           (recur (*' distance 2) (peek hare-path))))))
