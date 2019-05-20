(ns numthy.modular-arithmetic.abstract-algebra
  "Some basic functions exploring abstract algebra as it relates to modular arithmetic.
  They are merely meant to illustrate the expression of mathematic definitions in Clojure;
  they are not meant to be performant, and will run very slowly for inputs beyond small numbers."
  (:require [numthy.helpers :refer [coprime?]]
            [numthy.primes.is-prime :refer [prime?]]
            [numthy.primes.sieves :refer [coprimes-to]]))

;; http://mathonline.wikidot.com/algebraic-structures-fields-rings-and-groups

(defn Z
  "Given n, creates a set of ℤ_n, the integers mod n, with appropriate operations for
  mod-n arithmetic. Setting :coprimes to true will limit the set of elements to
  coprimes of n ∈ ℤ_n."
  [n & {:keys [coprimes] :or {coprimes nil}}]
  {:elements                (if coprimes
                              (set (filter #(coprime? % n) (range n)))
                              (set (range n)))
   :addition                (fn [a b] (mod (+ a b) n))
   :multiplication          (fn [a b] (mod (* a b) n))
   :additive-identity       0
   :multiplicative-identity 1})

(defn Zi
  "Given n, creates a set of ℤ_n[i], the Gaussian integers mod n, with appropriate
  operations for mod-n arithmetic over the Gaussian integers."
  [n]
  ;; (x+yi)(u+vi) = (xu-yv)+(xv-yu)i <-- https://www2.clarku.edu/faculty/djoyce/complex/mult.html
  ;; TODO: GCD/coprimality of Gaussian integers
  (let [rng (range n)]
    {:elements                (set (for [a rng b rng]
                                     {:re a :im b}))
     :addition                (fn [{re1 :re im1 :im} {re2 :re im2 :im}]
                                {:re (mod (+ re1 re2) n)
                                 :im (mod (+ im1 im2) n)})
     :multiplication          (fn [{re1 :re im1 :im} {re2 :re im2 :im}]
                                {:re (-> (* re1 re2) (- (* im1 im2)) (mod n))
                                 :im (-> (* re1 im2) (+ (* re2 im1)) (mod n))})
     :additive-identity       {:re 0 :im 0}
     :multiplicative-identity {:re 1 :im 0}}))

;; TODO: ZX --> ℤ_n[X], the polynomials mod n; dual numbers
;; https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields

(defn closed?
  "A set S is closed if a⋅b ∈ S ∀ a, b ∈ S."
  [m op]
  (let [S (:elements m) op (op m)]
    (every? true? (for [a S b S]
                    (contains? S (op a b))))))

(defn associative-set?
  "A set S is associative if (a⋅b)⋅c = a⋅(b⋅c) ∀ a, b, c ∈ S.
  (`-set` is appended to prevent overlap with `associative?` in clojure.core.)"
  [m op]
  (let [S (:elements m) op (op m)]
    (every? true? (for [a S b S c S]
                    (= (op (op a b) c)
                       (op a (op b c)))))))

(defn identity-element-naive
  "Is there an element e ∈ S s.t. a = a⋅e = e⋅a ∀ a ∈ S? If so, what is it?"
  [m op]
  (let [S (:elements m) op (op m)]
    (first (filter (fn [e]
                     (every? (fn [a] (= a (op a e) (op e a))) S))
                   S))))

(comment
 "In practice, we already know the identity element beforehand, so just do (:identity-element m)")

(defn identity-element
  ([m]
   (or (:additive-identity m) (:multiplicative-identity m)))
  ([m op]
   (condp = op
     :addition (:additive-identity m)
     :multiplication (:multiplicative-identity m)
     nil)))

(defn inverse
  "Given a, e ∈ S, is there an element a' ∈ S s.t a⋅a' = e, where is the identity element?
  If so, what is it?"
  [a m op]
  (when (contains? (:elements m) a)
    (let [e (identity-element m op)
          op (op m)]
      (first (filter (fn [a'] (= (op a a') e)) (:elements m))))))

(defn every-element-has-inverse?
  "Is there an element a' ∈ S s.t. a⋅a' = e, where e is the identity element, ∀ a ∈ S?"
  [m op]
  (every? (fn [a] (some? (inverse a m op))) (:elements m)))

(defn monoid?
  "A non-empty set S is a monoid if it is associative and has an identity element
  under some binary operation."
  [m op]
  (and (associative-set? m op)
       (some? (identity-element m op))))

(defn group?
  "A non-empty set S is a group under some binary operation if it is closed,
  associative, has an identity element, and every element has an inverse."
  [m op]
  (every? true? ((juxt closed?
                       associative-set?
                       (comp some? identity-element)
                       every-element-has-inverse?)
                 m op)))

(defn abelian?
  "An abelian group is a group that is commutative, i.e., a⋅b = b⋅a ∀ a, b ∈ S."
  [m op]
  (when (group? m op)
    (let [S (:elements m) op (op m)]
      (every? true? (for [a S b S c S]
                      (= (op a b) (op b a)))))))

(defn distributive?
  "A set S is distributive if a⋅(b+c) = (a⋅b)+(a⋅c) and (b+c)⋅a = (b⋅a)+(c⋅a) ∀ a, b, c ∈ R."
  [m]
  (let [S (:elements m) add-op (:addition m) mul-op (:multiplication m)]
    (every? true? (for [a S b S c S]
                    (and (= (mul-op a (add-op b c))
                            (add-op (mul-op a b) (mul-op a c)))
                         (= (mul-op (add-op b c) a)
                            (add-op (mul-op b a) (mul-op c a))))))))

(defn ring?
  "A set S is a ring if it is an abelian group under addition, a monoid under
  multiplication, and its multiplication is distributive over addition."
  [m]
  (and (abelian? m :addition)
       (monoid? m :multiplication)
       (distributive? m)))

(defn idempotents
  "The idempotent elements of a monoid, group, or ring are those a for which a⋅a=a.
  0 is always idempotent. 1 is idempotent under multiplication in ℤ. In ℤ_n, if
  n is composite, there are more idempotent elements than 0 and 1."
  [m op]
  (let [op (op m)]
    (set (filter (fn [a] (= (op a a) a)) (:elements m)))))

(defn integral-domain?
  "An integral domain is a ring R where a⋅b ≠ 0 ∀ a, b ≠ 0 ∈ R (0 denotes the zero
  element). That is, it is impossible to for two nonzero elements to `cancel out`
  through multiplication. For example, ℤ_7 is an integral domain, but ℤ_6 is not,
  because 2×3 = 0."
  [m]
  (when (ring? m)
    (let [R        (:elements m)
          nonzero? #(not= % (:additive-identity m))
          mul-op   (:multiplication m)]
      (every? nonzero? (for [a R b R :when (and (nonzero? a) (nonzero? b))]
                         (mul-op a b))))))

(defn finite-field?
  "A set S is a field if it is commutative and associative under addition and multiplication,
  it has an additive and multiplicative identity, each of its elements has an additive inverse
  (equivalent to subtraction), each of its nonzero elements has a multiplicative inverse
  (equivalent to division), and its multiplication is distributive over its addition. A finite
  field is simply a field with finitely many elements."
  [m]
  (let [nonzero-m (update m :elements disj (:additive-identity m))]
    (and (abelian? m :addition)
         (abelian? nonzero-m :multiplication)
         (distributive? m))))
