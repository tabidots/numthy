(ns numthy.factorization.core
  (:require [clojure.core.async :as a
             :refer [>! <! >!! <!! go go-loop chan sliding-buffer alts! alts!! close! thread timeout]]
            [numthy.factorization.mpqs :refer [mpqs]]
            [numthy.factorization.pollard :refer [rho brent fast-brent]]
            [numthy.factorization.squfof :refer [squfof]]
            [numthy.helpers :refer [count-bits]]
            [numthy.perfect-powers :refer [perfect-power?]]
            [numthy.primes.is-prime :refer [quick-prime?]]
            [numthy.primes.sieves :refer [primes-to]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

(def tdiv-primes (primes-to 100000))

(defn- divide-out-factors
  "Divides n over each x in coll, returning a hash-map with the remainder (an integer)
  and a hash-map of {x_i e_i}, where e is the maximal power of x that divides n."
  [n coll]
  (let [n' (volatile! n) m (volatile! {})]
    (doseq [p coll]
      (while (zero? (mod @n' p))
        (vswap! n' quot p)
        (vswap! m update p (fnil inc 0))))
    {:remainder @n' :factors @m}))

(defn- divide-out-small-factors
  "Divides out prime factors < 100000 from n, returning a map with the remainder and
  factorization."
  [n]
  (do (println "1. Dividing out small factors...")
    (->> (take-while #(<= % n) tdiv-primes)
         (filter #(zero? (mod n %)))
         (divide-out-factors n))))

(defn- divide-out-medium-factors
  "If n is still not completely factored after dividing out factors < 100000,
  uses the SQUFOF and Pollard-Brent algorithms to find larger factors. SQUFOF is
  only engaged for remainders up to 64 bits. Pollard-Brent is limited here to about
  14 iterations to keep the running time down."
  [{:keys [remainder factors] :as m}]
  (if (= 1 remainder) m
    (do
      (println "2. Dividing out medium-sized factors with SQUFOF and fast Pollard-Brent...")
      (loop [n remainder res factors]
        (if (quick-prime? n)
          {:remainder 1 :factors (merge-with + res {n 1})}
                ;; 3 dice-rolls of Brent may find 0, 1 or 2 distinct factors
          (let [brents (filter quick-prime? (take 3 (repeatedly #(fast-brent n))))
                fs     (remove nil? (conj brents (when (<= (count-bits n) 64)
                                                   (squfof n))))]
            (if (empty? fs)
              {:remainder n :factors res} ;; Factoring as complete as possible at this stage
              (let [{n' :remainder m :factors} (divide-out-factors n (distinct fs))]
                (recur n' (merge res m))))))))))

(defn- divide-out-large-factors
  "If n is still not completely factored after dividing out medium-sized factors,
  uses Pollard-Brent (up to about 21 iterations) and MPQS to find larger factors."
  [{:keys [remainder factors] :as m}]
  (if (= 1 remainder) m
    (let [result-chan  (chan (sliding-buffer 1))
          timeout-chan (timeout 60000)]
      (println "3. Looking for large factors with Pollard-Brent and MPQS...")
      (loop [n remainder res factors]
        (if (quick-prime? n)
          (do (println "Finished!")
            (close! result-chan)
            {:remainder 1 :factors (merge-with + res {n 1})})
          ;; This approach is really ugly, but it's the only way to actually definitively interrupt
          ;; the processes in the middle of their execution, to keep them from hogging the CPU forever
          ;; https://stackoverflow.com/a/42700636/4210855
          (let [mpqs-thread  (Thread. (fn [] (try
                                               (while true
                                                 (some->> (mpqs n) (>!! result-chan)))
                                               (catch InterruptedException e
                                                 (println (.getMessage e))))))
                brent-thread (Thread. (fn [] (try
                                               (while true
                                                 (some->> (brent n) (>!! result-chan)))
                                               (catch InterruptedException e
                                                 (println (.getMessage e))))))]
            (.start mpqs-thread)
            (.start brent-thread)
            (let [[p c] (alts!! [result-chan timeout-chan])]
              (.interrupt mpqs-thread)
              (.interrupt brent-thread)
              (if (= c timeout-chan)
                (do (println "Timed out. This is as far as we could get:")
                  {:remainder n :factors res})
                (do (println "Found a factor! Continuing...")
                  (recur (quot n p) (merge-with + res {p 1})))))))))))

(defn factorize
  "Returns the prime factorization of n as a hash-map of the form {p_i e_i}."
  [n]
  (when (> n 1)
    (if (quick-prime? n) {n 1}
      (if-some [pp (perfect-power? n)] pp
        (->> (divide-out-small-factors n)
             (divide-out-medium-factors)
             (divide-out-large-factors)
             :factors
             (into (sorted-map)))))))

(defn distinct-prime-factors
  "Returns the distinct prime factors of an integer n in a vector. e.g., 1300 = [2 5 13]"
  [n]
  (some->> (factorize n) keys vec))

(defn prime-factors-with-multiplicity
  "Returns a vector of the prime factors of an integer n with multiplicity.
  e.g., 168 = [2 2 2 2 3 7]"
  [n]
  (some->> (factorize n)
           (reduce-kv (fn [r b e]
                        (into r (repeat e b)))
                      [])))

(defn prime-signature
  "Returns the non-zero exponents of the prime factorization of n, sorted in a vector."
  [n]
  (some->> (factorize n) vals sort vec))

(defn phi
  "Uses the product rule to compute Euler's totient function of n, which the number of
  integers less than n that are coprime to it."
  [n]
  (->> (distinct-prime-factors n)
       (map #(- 1 (/ 1 %)))
       (reduce * n)))
