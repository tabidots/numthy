(ns numthy.primes.sieves
  (:require [numthy.helpers :refer [isqrt]]))

;; TODO: Test performance with type hints and specific fns (aset-boolean, etc)

(defn primes-to
  "Finds all prime numbers less than n and returns them sorted in a vector."
  ;; Source: https://stackoverflow.com/a/48137788/4210855
  [n]
  (if (< n 2) []
    (let [^booleans sieve (boolean-array n false)
          s (-> n Math/sqrt Math/floor int)]
      (loop [p 2]
        (if (> p s)
          (into []
            (remove #(aget sieve %))
            (range 2 n))
          (do
            (when-not (aget sieve p)
              (loop [i (* 2 p)]
                (when (< i n)
                  (aset sieve i true)
                  (recur (+ i p)))))
            (recur (inc p))))))))

(defn spf-sieve
  "A modified Sieve of Eratosthenes, using a JVM int-array, that records the
  smallest prime divisor of all numbers < n, rather than just prime?=true/false.
  Does not drop the first two values (for 0 and 1)."
  ;; Adapted from https://stackoverflow.com/a/48137788/4210855
  [n]
  (let [^ints sieve (int-array n (range n))
        upper-bound (isqrt n)]
    (loop [p 2] ;; p's are the bases
      (if (> p upper-bound) sieve
        (do
          (when (= p (aget sieve p))
            (loop [i (* 2 p)] ;; i's are the multiples of p
              (when (< i n)
                (aset sieve i (min (aget sieve i) p)) ;; Don't overwrite with larger p's
                (recur (+ i p)))))
          (recur (inc p)))))))

(defn squarefrees-to
  "Uses a modified Sieve of Eratosthenes to find all square-free numbers less
  than n and returns them sorted in a vector."
  ;; Source: https://stackoverflow.com/a/48137788/4210855
  [n]
  (if (< n 2) [1]
    (let [^booleans sieve (boolean-array n false)
          s (-> n Math/sqrt Math/floor int)]
      (loop [p 2]
        (if (> p s)
          (into []
                (remove #(aget sieve %))
                (range 1 n))
          (do
            (when-not (aget sieve p)
              (let [sq (*' p p)]
                (loop [i sq]
                  (when (< i n)
                    (aset-boolean sieve i true)
                    (recur (+ i sq))))))
            (recur (inc p))))))))
