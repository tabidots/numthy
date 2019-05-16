(ns numthy.factorization.benchmark
  "Functions to generate test suites for factorization algorithms."
  (:require [clj-async-profiler.core :as prof]
            [clojure.core.reducers :as r]
            [clojure.main :refer [demunge]]
            [clojure.string :as s]
            [criterium.core :as c]
            [criterium.stats :refer [mean median]]
            [numthy.factorization.pollard :refer [rho brent]]
            [numthy.factorization.cfrac :refer [cfrac]]
            [numthy.factorization.dixon :refer [dixon]]
            [numthy.factorization.mpqs :refer [mpqs]]
            [numthy.factorization.squfof :refer [squfof]]
            [numthy.helpers :refer [divisible? count-bits count-digits]]
            [numthy.primes.generate :refer [rand-prime]]
            [numthy.perfect-powers :refer [perfect-power?]])
  (:import java.util.concurrent.ThreadLocalRandom))

(prof/serve-files 8080)

(def log2-10
  "The binary logarithm of 10, which, when multiplied by a number d, gives
  the maximum number of bits needed for a d-digit number."
  ;; https://www.exploringbinary.com/number-of-bits-in-a-decimal-integer/
  (with-precision 75 (/ (BigDecimal. (Math/log 10)) (BigDecimal. (Math/log 2)))))

(defn rand-semiprime
  "Generates a random n with the specified number of bits, s.t. n = pq with p and q prime
  and both somewhat close to √n - i.e., a composite that is theoretically more difficult to
  factor than general composites with the same number of bits."
  [bits1 bits2]
  (let [p (rand-prime bits1)
        q (->> #(rand-prime bits2)
               repeatedly
               (remove #{p})
               first)]
    (*' p q)))

(defn n-digit-rand-semiprime
  "Generates a random n with the specified number of digits, s.t. n = pq with p and q prime
  and both somewhat close to √n - i.e., a composite that is theoretically more difficult to
  factor than general composites with the same number of digits."
  [digits]
  (let [max-bits  (inc (int (*' digits log2-10)))
        half-bits (quot max-bits 2)]
    (rand-semiprime half-bits half-bits)))

(defn rand-easy-semiprime
  "Generates a random integer with the specified number of bits, s.t. n = pq with p and q prime and p < ∛n."
  [bits]
  (let [p-bits (dec (quot bits 3)) ;; dec for wiggle room
        q-bits (- bits p-bits)]
    (rand-semiprime p-bits q-bits)))

(defn n-digit-rand-easy-semiprime
  "Generates a random integer with the specified number of digits, s.t. n = pq with p and q prime and p < ∛n."
  [digits]
  (let [max-bits (inc (int (*' digits log2-10)))
        p-bits   (dec (quot max-bits 3)) ;; dec for wiggle room
        q-bits   (- max-bits p-bits)]
    (rand-semiprime p-bits q-bits)))

(defn rand-good-composite
  [bits]
  (let [src (ThreadLocalRandom/current)]
    (->> #(BigInteger. ^long bits src)
         repeatedly
         (remove (fn [n]
                   (or (< n 1000)
                       (.isProbablePrime ^BigInteger n 5)
                       (some (fn [x] (divisible? n x)) (range 2 (Math/log10 n)))
                       (perfect-power? n))))
         first)))

(defn generate-distinct
  [how-many what]
  (r/reduce (fn [res n]
              (if (>= (count res) how-many) (reduced res)
                (conj res n)))
            #{} (repeatedly what)))

(defn fn-name
  [f]
  (re-find #"(?<=/)\S*(?=@)" (demunge (str f))))

(defn convert-time
  [millis]
  (cond
    (>= millis 60000) (let [[mins secs] ((juxt quot rem) millis 60000)]
                        (str (int mins) " min " (convert-time secs)))
    (>= millis 1000)  (let [[secs mils] ((juxt quot rem) millis 1000)]
                        (str (int secs) " sec " (int mils) " ms"))
    :else             (str millis " ms")))

(defn run-test
  [sample-size generator & algos]
  (let [samples (generate-distinct sample-size generator)]
    (println "Samples generated.")
    (println (format "Mean bit size of samples, rounded down: ~%d bits." (int (mean (map count-bits samples)))))
    (doseq [algo algos]
      (println "")
      (println (format "Factoring with %s..." (fn-name algo)))
      (let [times* (transient {})
            fail   (volatile! 0)]
        (doseq [n samples]
          (let [f    (volatile! nil)

                s    (with-out-str (time (vreset! f (algo n))))
                stat (read-string (re-find #"\d[\S]*" s))]
            (assoc! times* n stat)
            (when (nil? @f)
              (vswap! fail inc))))
        (let [times      (persistent! times*)
              failures   @fail
              worst-case (apply max-key val times)]
          (println (format "Median time: %s" (convert-time (first (median (vals times))))))
          (println (format "  Mean time: %s" (convert-time (mean (vals times)))))
          (println (format " Worst case: %s (n = %d)" (convert-time (val worst-case))
                                                      (biginteger (key worst-case))))
          (println (format "   Failures: %d (Success rate: %f%%)" failures
                                                                  (- 100 (* 100.0 (/ failures sample-size))))))))))

(defn test-semiprimes [sample-size digits & algos]
  (println (format "Test: Factoring %d %d-digit semiprimes." sample-size digits))
  (println (format "Generating %d random samples..." sample-size))
  (apply run-test sample-size #(n-digit-rand-semiprime digits) algos))

(defn test-composites [sample-size bits & algos]
  (println (format "Test: Factoring %d %d-bit composites." sample-size bits))
  (println (format "Generating %d random samples..." sample-size))
  (apply run-test sample-size #(rand-good-composite bits) algos))
