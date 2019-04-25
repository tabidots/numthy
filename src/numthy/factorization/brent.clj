(ns numthy.factorization.brent
  (:require [clojure.math.numeric-tower :refer [gcd expt abs]]
            [numthy.helpers :refer [rand-num]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.modular-arithmetic.utils :refer [mod-pow mod-mul]]
            [clojure.core.reducers :as r]))

(defn pollard-brent
  "Uses the Pollard-Brent algorithm to find one prime factor of a large composite integer n."
  [n]
  (if (even? n) 2
    (let [c (rand-num 1 (dec n))                       ;; c = constant of the polynomial
          m (rand-num (Math/log n) (if (< 2000) (dec n) ; m = chunk size
                                     (expt n 1/4)))
          f (fn [x] (-> (+ c (mod-pow x 2 n)) (mod n)))
          x (rand-num 1 (dec n))]
      (loop [r 1 q 1 x x]                              ;; x  = Tortoise value
        (let [y   (f x)                                ;; y  = Backup value if algo fails
              ys  (into [] (r/take r (iterate f y)))   ;; ys = Sequence of hare values
              qs  (->> (r/map #(abs ^BigInteger (- x %)) ys)
                       (into [])
                       (reductions (fn [res z]         ;; qs = Cumulative products of tort/hare distances
                                     (mod-mul res z n)) ;      all the way from the beginning
                                   q)
                       vec)
              q*  (peek qs)
              g   (or (->> (take-nth m qs)              ;; We want to compute g = gcd(q, n) in chunks, at
                           (r/map #(gcd % n))           ;; (1) multiples of the step size, hence take-nth, or
                           (r/filter #(> % 1))          ;; (2) a value of r, hence cons the last value.
                           (r/foldcat)                  ;;
                           first)                       ;; If all of those values of g are 1, leave g = 1 so the
                      1)]                               ;; loop continues. Otherwise, we found a factor
          (cond
            (< 1 g n) g                                 ;; Success!
            (= g n)   (->> ys                           ;; Second attempt, iterating in steps of 1 instead of m
                           (r/map #(gcd ^BigInteger (abs ^BigInteger (-' x %)) n))
                           (r/filter #(< 1 % n))
                           (r/foldcat)
                           first)
            ;; If r gets too high, we may be dealing with a semi prime, or otherwise very large prime factors.
            ;; (> r ____)  (some-other-factorization-method)
            :else     (recur (*' 2 r) q* (peek ys))))))))

(defn brent-factorize
  "Uses the Pollard-Brent algorithm repeatedly to find all prime factors of a large
  integer n, along with their powers, in a hash-map of the form {p_0 e_0, ⋯, p_i e_i}."
  [n]
  (when (and (integer? n) (> n 1))
    (loop [n n res (sorted-map)]
      (if (.isProbablePrime (biginteger n) 5)
        (merge-with + res {n 1})
        (let [p (first (filter #(.isProbablePrime (biginteger %) 5)
                               (repeatedly #(pollard-brent n))))]
          (recur (/ n p) (merge-with + res {p 1})))))))

;; Some test values

;; 209348923723492392323490230948833189789374897823432748973
;; 5687450887724258126853054080419196271544306051 is a composite factor of that, takes a long time to factorize

;; c 1384618644036103387726606019273596888099126033
;; m 163544
;; x 117210326061389040030818806148762920404935678

;; 209348923723492392323490230948833189789374897823432748973
;; 4493699349619351139582065363011367885289N is a composite factor of that, takes a long time to factorize
