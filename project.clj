(defproject numthy "0.1.0-SNAPSHOT"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :jvm-opts ^:replace ["-Djdk.attach.allowAttachSelf" "-server"]
  :exclusions [[org.jcuda/jcuda-natives :classifier "apple-x86_64"]
               [org.jcuda/jcublas-natives :classifier "apple-x86_64"]]
  :dev-dependencies []
  :dependencies [[org.clojure/clojure "1.10.0"]
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [org.clojure/math.combinatorics "0.1.4"]
                 [net.mikera/vectorz-clj "0.48.0"]
                 [net.mikera/core.matrix "0.62.0"]
                 [org.clojure/core.async "0.4.490"]
                 [org.flatland/useful "0.11.6"]
                 [criterium "0.4.4"]
                 [com.clojure-goes-fast/clj-async-profiler "0.3.1"]])
