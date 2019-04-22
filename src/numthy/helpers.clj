(ns numthy.helpers
  (:require [clojure.string :as s]
            [clojure.math.numeric-tower :refer [gcd expt]]
            [clojure.math.combinatorics :refer [combinations]])
  (:import java.util.concurrent.ThreadLocalRandom))

;; GENERAL UTILITIES

(defn divisible? [dividend divisor]
  (zero? (mod dividend divisor)))

(defn coprime? [a b]
  (= 1 (gcd a b)))

(defn pairwise-coprime?
  "True if every distinct pair of (a, b) in xs is coprime."
  [xs]
  (every? (partial apply coprime?) (combinations xs 2)))

(defn pfilter
  "Like filter, except uses multiple cores. Useful for huge collections."
  [pred coll]
  (->> (pmap #(when (pred %) %) coll)
       (remove nil?)))

(defn rand-num
  "Random number generator of arbitrary-precision integers."
  ([max]
   (rand-num 0 max))
  ([min max]
   (let [src (ThreadLocalRandom/current)]
     (cond
       ;; Clojure translation of https://stackoverflow.com/a/2290089
       (> max Long/MAX_VALUE)    (->> #(BigInteger. (.bitLength (bigint max)) src)
                                      repeatedly
                                      (filter #(<= min % max))
                                      first)
       ;; For longs instead of ints, use ThreadLocalRandom.current().nextLong(start, end)
       ;; Clojure High Performance Programming, p. 75
       (> max Integer/MAX_VALUE) (.nextLong src min max)
       :else                     (.nextInt src min max)))))

(defn isqrt
  "floor(√n). When incremented, provides an upper bound for factorization."
  ;; Java interop is super fast but not accurate for n > 1E24 (approx) due to
  ;; floating-point rounding. Uses a slightly slower but pinpoint-precise method for n > 1E24.
  [n]
  (if (< n 1E24)
    (-> (Math/sqrt n) bigint)
    ;; https://cs.stackexchange.com/a/30383
    (let [half-bit-length (quot (.bitLength (bigint n)) 2)]
      (loop [a (expt 2 half-bit-length)
             b a
             c (*' a a)]
        (cond
          (zero? b) a
          (> c n)   (recur (-' a b) (quot b 2) (+ c (*' -2 a b) (*' b b)))
          :else     (recur (+' a b) (quot b 2) (+ c (*' 2 a b) (*' b b))))))))

(defn factorial [n]
  (reduce * (take n (iterate (partial inc) (bigint 1)))))

;;;;;;;;;

(comment
 (defn phi
   "Euler's totient function, optimized using the product rule."
   [n]
   ;; https://en.wikipedia.org/wiki/Euler%27s_totient_function#Example
   (let [prime-factors (filter prime? (factors n))
         diffs         (map #(- 1 (/ 1 %)) prime-factors)]
     (reduce * n diffs)))

 (defn prime-factorization [n]
   "Returns the prime factorization of an integer, e.g., 168 -> [2 2 2 3 7]."
   (when (and (integer? n) (pos? n))
     (loop [n       n
            primes  (if (> n 1000000)
                      (filter prime? (factors n))
                      (filter prime? (range)))
            factors []]
       (let [k (first primes)]
         (cond
           (= 1 n)          factors
           (divisible? n k) (recur (/ n k) primes (conj factors k))
           :else            (recur n (rest primes) factors)))))))

; https://en.wikipedia.org/wiki/Exponentiation_by_squaring
(defn exp-by-sq [x n]
  (cond
    (neg? n)   (exp-by-sq (/ 1 x) (- n))
    (zero? n)  1
    (= 1 n)    x
    (even? n)  (exp-by-sq (* x x) (/ n 2))
    (odd? n)   (exp-by-sq (* x x) (-> (- n 1) (/ 2)))))

;; SEQUENCES

(defn fibonacci
  ([]
   (->> (iterate (fn [[f0 f1]] [f1 (+ f0 f1)]) [0N 1N])
        (map peek)))
  ([n]
   (take n (fibonacci))))
