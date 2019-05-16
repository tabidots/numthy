(ns numthy.factorization.mpqs
  "Multiple-Polynomial Quadratic Sieve."
  (:require [clojure.core.matrix :refer [order]]
            [clojure.core.reducers :as r]
            [clojure.math.numeric-tower :refer [expt abs gcd]]
            [flatland.useful.map :refer [map-to]]
            [numthy.factorization.squares-utils :refer [find-winning-subsets smoothness-bound]]
            [numthy.helpers :refer [isqrt]]
            [numthy.modular-arithmetic.utils :refer [mod-inverse mod-mul mod-pow]]
            [numthy.modular-arithmetic.quadratic-residue :refer [msqrt]]
            [numthy.polynomials.core :refer [poly-eval]]
            [numthy.primes.is-prime :refer [quick-prime?]]
            [numthy.primes.sieves :refer [primes-to]]))

(set! *unchecked-math* true)
(set! *warn-on-reflection* true)

;; https://www.rose-hulman.edu/~bryan/lottamath/factor.pdf

;; This isn't currently being used, but keeping it here just in case
(defn qs-smoothness-bound
  "Floor of e ^ ((√2/4) * √(ln n ln ln n)), as per Landquist's QS implementation."
  [n]
  (let [sqrt2-over-4 (/ (Math/sqrt 2) 4)
        lognloglogn  (*' (Math/log n) (Math/log (Math/log n)))]
    (biginteger (Math/exp (* sqrt2-over-4 (Math/sqrt lognloglogn))))))

;; UTILITY FUNCTIONS

(defn fast-legendre
  "Short helper function extending the Legendre symbol to the case p = 2, and also
  without the primality test because the only numbers fed into this will be primes.
  Also omits the modulo wraparound, so returns 1 if a is a quadratic residue mod p,
  0 if a divides p, and (p-1) if a is not a quadratic residue mod p."
  [a p]
  (if (= p 2) (if (odd? a) 1 0)
    (mod-pow a (-> (- p 1) (/ 2)) p)))

(defn next-prime
  "Quickly generates the first prime greater than n."
  [n]
  (if (< n 2) 2
    (loop [p (inc n)]
      (if (quick-prime? p) p
        (recur (inc p))))))

;; MAJOR COMPONENTS

(defn make-factor-base
  "Given an integer n to be factored, calculates a smoothness bound B and returns a
  vector of at least 3 primes of which n is a quadratic residue. B is the upper limit for
  the highest prime, unless fewer than 3 eligible primes are found by that point,
  in which case more primes will be tested.
  NOTE: If a non-trivial factor of n is found before the factor base is populated,
  return early."
  [n]
  (let [B           (* 5 (expt (Math/log10 n) 2))
        wiggle-room 1.2] ;; Multiply B by some number to ensure we don't run out of eligible primes
    (loop [base   (transient [])
           primes (primes-to (biginteger (* wiggle-room B)))]
      (if-let [p (first primes)]
        (let [ls (fast-legendre n p)]
          (cond
            (and (>= (count base) 3)
                 (>= p B))  {:n n :primes (persistent! base)}
            (zero? ls)      {:factor p}  ;; Trial div found a factor
            (not= 1 ls)     (recur base (rest primes))
            :else           (recur (conj! base p) (rest primes))))
        (println "Ran out of primes!"))))) ;; For debugging: Increase wiggle room if this happens too often

(defn initialize-mpqs
  "Initializes a Multiple-Polynomial Quadratic Sieve to factor n, given a factor base."
  [base]
  (let [{n      :n
         primes :primes} base
        sqrt2n           (isqrt (+' n n))
        sieve-radius     (* 60 (count primes))               ;; corresponds to "M"
        max-gx           (*' sieve-radius (/ sqrt2n 2)) ;; maximum possible g(x)
        base-threshold   (* (Math/log10 max-gx) 0.735)
        quadratic-roots  (map-to #(or (msqrt n %)
                                      (if (odd? n) 1 0)) ;; safe for mod 2
                                 (rest primes))
        ;; Optimization #1: "Small Prime" variant
        ;; Define a minimum value s.t. lesser primes are not explicitly sieved,
        ;; because their contribution is not worth the time spent in computation,
        ;; but rather the log threshold is adjusted downwards a bit to compensate
        min-sig-prime    (int (* base-threshold 3))
        fudge-factor     (->> primes
                              (filter #(< % min-sig-prime))
                              (pmap #(Math/log10 %))
                              (reduce +)
                              (* 1/4))
        log-threshold    (- base-threshold fudge-factor)]
    (merge base {:sqrt2n sqrt2n
                 ;; For generating polynomials
                 :seed (isqrt (/ sqrt2n sieve-radius)) ;; corresponds to q in Contini
                 ;; For sieving
                 :M sieve-radius :log-threshold log-threshold
                 :min-sig-prime min-sig-prime :quadratic-roots quadratic-roots
                 ;; For collecting relations
                 :partials {} :smooths (sorted-map)
                 ;; For generating the matrix
                 :used-primes (sorted-set)})))

(defn next-polynomial
  "Generates a polynomial that meets the requirements of the MPQS based on a generated
  or iterated seed value."
  [state]
  (let [{n :n seed :seed roots :quadratic-roots} state]
    (loop [sqrt-A (next-prime seed)]
      (condp = (fast-legendre n sqrt-A)
        0 {assoc state :factor sqrt-A} ;; Stumbled upon a factor
        1 (let [A     (*' sqrt-A sqrt-A)
                b     (or (msqrt n sqrt-A)
                          (if (odd? n) 1 0)) ;; safe for mod 2
                x1    (or (mod-inverse (+ b b) sqrt-A)
                          0) ;; safe for mod 2
                x2    (- n (* b b))
                B     (-> (* x1 x2) (+ b) (mod A))
                C     (-> (*' B B) (- n) (quot A))
                solns (reduce-kv (fn [res p t]
                                   (if-some [A-1 (mod-inverse A p)]
                                     (assoc res p
                                            [(-> (- t B) (* A-1) (mod p))
                                             (-> (- p t B) (* A-1) (mod p))])
                                     (assoc res p [0])))
                                 {} roots)]
            (merge state {:seed sqrt-A
                          ;; The full polynomial is (Ax + B)^2. But we've purposely chosen A to be a square, so
                          ;; we can factor out A and just consider Ax^2 + 2Bx + C -> {2 A, 1 2B, 0 C}
                          :polynomial {2 A, 1 (+' B B), 0 C}
                          :sqrt-polynomial {1 A, 0 B} ;; We still need the square root of the full polynomial
                          :quadratic-solutions solns}))
        (recur (next-prime sqrt-A))))))

(defn sieve-candidates
  "Performs a quadratic sieve over the interval [-M, M] by adding the log10 of each prime
  to values of x in the interval that are divisible by g(x) and filtering the ones whose final
  log sum is greater than the threshold value. Then saves the x and g(x) values of each
  candidate in a map, and puts them all in a vector."
  [state]
  (if (some? (:factor state)) state
    (let [{min-sig-prime   :min-sig-prime
           primes          :primes
           solutions       :quadratic-solutions
           pnml            :polynomial
           threshold       :log-threshold
           M               :M}  state
          ^doubles log-sums     (double-array (inc (* 2 M)) 0.0) ;; 2M+1 = -M to M, incl. zero
          candidates            (transient [])]
      ;; Easy breezy approximate "trial division" through accumulation of logarithms of significant primes at
      ;; regular points throughout the interval as indicated by the modular square roots of those primes
      (doseq [p (filter #(>= % min-sig-prime) primes)] ;; See note about "Small Prime" variant above
        (let [logp (Math/log10 p)]
          (doseq [s (get solutions p)]
            ;; Array is meant to reflect values from -M to M, but is zero-indexed
            ;; Need to find the index of the first value of x after -M s.t. x ≡ s mod p
            ;; Looping or using `first` to find starting points is slow b/c repeated mod calculations
            ;; This is an ugly but fast & efficient way to find starting point with only one mod calculation
            ;; start := (-M mod p) - s + p
            (loop [i (- p (- (mod (- M) p) s))]
              (if (Thread/interrupted)
                (throw (InterruptedException. "Interrupting MPQS..."))
                (when (<= i (* 2 M))
                  (aset ^doubles log-sums i (+ (aget ^doubles log-sums i) logp))
                  (recur (+ i p))))))))
      (dotimes [i (inc (* 2 M))] ;; This can be written functionally as keep-indexed, but is much slower
        (let [log-gx (aget ^doubles log-sums i)]
          (when (> log-gx threshold)
            (let [x (- i M) gx (poly-eval pnml x)]
              (conj! candidates {:gx gx :x x})))))
      (assoc state :candidates (persistent! candidates)))))

(declare check-smooth)

(defn add-relations
  "Tests each candidate (value of g(x)) and adds it to the map `smooths` if it is
  smooth over the factor base, `partials` if it is smooth except for one factor, and
  merges any matching partials into a smooth."
  [state]
  (if (some? (:factor state)) state
    (let [{candidates :candidates primes :primes} state]
      (-> (reduce (fn [state' candidate]
                    (if (Thread/interrupted)
                      (throw (InterruptedException. "Interrupting MPQS..."))
                      (->> (check-smooth state' candidate)
                           (merge-with into state'))))
                  state candidates)
          (dissoc :candidates)))))

(defn common-factors
  "Finds common factors between two exponent maps."
  [m1 m2]
  (reduce-kv (fn [r k v]
               (if (= 1 v (get m2 k 0))
                 (conj r k)
                 r))
             [] m1))

(defn check-smooth
  "Checks one candidate by trial division to see if it is smooth over the factor base.
  Saves relevant information and adds it to `smooths` if it is smooth, adds it to `partials`
  if it is smooth but for one factor, or combines it with a matching partial if there is one."
  [state candidate]
  (let [{sqrt-A     :seed
         lhs        :sqrt-polynomial
         partials   :partials} state ;; partials for optimiziation #2: "Double Large Prime" variant.
        gx              (:gx candidate)
        square-factors  (volatile! [])
        used-primes     (volatile! (if (neg? gx) [-1] []))
        Ax+B            (poly-eval lhs (:x candidate))]
    (loop [gx'          (abs gx)
           primes       (:primes state)
           exponent-map (if (neg? gx) {-1 1} {})]
      (if (= 1 gx') ;; Candidate is smooth
        {:smooths     {Ax+B {:exponent-map   exponent-map
                             :square-factors @square-factors
                             :rhs-cofactor   sqrt-A}}
         :used-primes @used-primes}
        (if-some [p (first primes)]
          (if (zero? (mod gx' p))
            (do
              (when (= 1 (get exponent-map p)) ;; p is a square factor of gx'
                (vswap! square-factors conj p))
              (vswap! used-primes conj p)
              (recur (/ gx' p) primes (update exponent-map p (fnil bit-xor 0) 1)))
            (recur gx' (rest primes) exponent-map))
          (if-some [partial (get partials gx')]
            ;; If no primes left, and another candidate shares the same remaining value,
            ;; combine the partials to make a smooth.
            {:smooths {(*' Ax+B (:lhs partial))
                       {:exponent-map   (merge-with bit-xor exponent-map (:exponent-map partial))
                        :square-factors (-> (concat @square-factors
                                                    (:square-factors partial)
                                                    (common-factors exponent-map (:exponent-map partial)))
                                            (conj gx'))
                        :rhs-cofactor   (*' sqrt-A (:rhs-cofactor partial))}}
             :used-primes @used-primes}
            ;; Otherwise, save gx' as a partial key
            {:partials    {gx' {:exponent-map   exponent-map
                                :lhs            Ax+B
                                :square-factors @square-factors
                                :rhs-cofactor   sqrt-A}}
             :used-primes @used-primes}))))))

(defn create-matrix
  "Creates a GF(2) matrix from a sequence of exponent maps, based on the primes that
  have been used. These may be different than the original factor base."
  [state]
  (if (some? (:factor state)) state
    (let [{used-primes :used-primes smooths :smooths} state]
      ;; Pass along the state if not enough relations, to force another iteration
      (if (<= (count smooths) (count used-primes)) state
        (let [blank-map (into (sorted-map) (zipmap used-primes (repeat 0)))]
          (->> (vals smooths)
               (mapv (comp vec
                           vals
                           #(select-keys % used-primes)
                           #(merge blank-map %)
                           :exponent-map))
               (assoc state :matrix)))))))

(defn no-factors?
  "Returns true if no factor has been found yet."
  [state]
  (nil? (:factor state)))

(defn make-y-from
  "Constructs the big side of a congruence of squares from information about
  factors of the candidates belonging to a subset of rows."
  [ms]
  (let [sqfs      (r/fold *' (r/mapcat :square-factors ms))
        cofacts   (r/fold *' (r/map :rhs-cofactor ms))
        base-prod (->> (pmap :exponent-map ms)
                       (apply merge-with +)
                       (r/fold (r/monoid *' (constantly 1))
                               (fn ([res b e]
                                    (*' res (expt b (quot e 2)))))))]
    (*' base-prod cofacts sqfs)))

(defn make-x-from
  "Constructs the small side of a congruence of squares from values taken from
  a subset of rows."
  [ks]
  (r/fold *' ks))

(defn find-factors
  "Tests subsets of a matrix that sum to the zero vector in order to find a
  congruence of squares x and y s.t. gcd(x - y) mod n is a non-trivial factor of n."
  [state]
  (if (or (some? (:factor state)) (nil? (:matrix state))) state
    (if (Thread/interrupted)
      (throw (InterruptedException. "Interrupting MPQS..."))
      (let [{n :n matrix :matrix smooths :smooths} state]
        (->> (for [subset (find-winning-subsets matrix)]
               (let [[x y] (->> (order (vec smooths) subset)
                                ((juxt (comp make-x-from keys) (comp make-y-from vals))))]
                 (gcd (- x y) n)))
             (filter #(< 1 % n))
             (first)
             (assoc state :factor))))))

(defn mpqs
  "Finds one prime factor of n using a Multiple-Polynomial Quadratic Sieve."
  [n]
  (let [base (make-factor-base n)]
    (if-some [factor (:factor base)] factor
      (->> (initialize-mpqs base)
           (iterate (comp find-factors
                          create-matrix
                          add-relations
                          sieve-candidates
                          next-polynomial))
           (drop-while no-factors?)
           first :factor))))
