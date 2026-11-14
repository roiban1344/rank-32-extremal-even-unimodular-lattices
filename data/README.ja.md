# data/

[English](README.md) | 日本語

`classified-lattices.jsonl` は 1 行に 1 個の JSON オブジェクトを持ち、Aut L ≠ {±1} である階数 32 の極値偶ユニモジュラー格子 L
143329 個のそれぞれに 1 行が対応する（定理 B）。

| フィールド | 意味 |
|---|---|
| `label` | `32.<\|Aut L\|>.<i>`。データベース [D] のラベルである。番号の付け方は論文の §6.8 で定める |
| `aut` | Aut L の位数 |
| `part` | 定理 B の部分: `general`（S_gen）、`c2xc2`（S_(C_2^2)）、`cyclic-quaternion`（S_CQ） |
| `group.small_group_id` | 同定されている場合は、GAP の Small Groups Library での `[order, number]`。そうでなければ `null` |
| `group.structure` | Aut L の構造の記述（記録がある場合）。S_gen の格子の大半（12546 個のうち 12466 個）では `null` であり、それらの群は `small_group_id` だけで（12447 個）、または主列で（9 個）同定されているか、まったく同定されていない（\|Aut L\| > 10^7 である 10 個の格子。下を見よ） |
| `group.chief_series_orders` | 主列の位数。この方法だけで同定された群について記す |
| `group.presentation_checked` | Aut L の表示をその生成元の上で検査し、位数を確かめたこと（S_gen の 12536 個の格子で `true`）。S_(C_2^2) と S_CQ（群は第 1 の論文の §§6.4–6.5 で決定される）、および \|Aut L\| > 10^7 である S_gen の 10 個の格子（表示は第 2 の論文で証明され、`largest-groups.jsonl` に収めてある）では `false` |
| `entry_types` | L の項目の種類（§6.1）を文字列で記す: 対合の項目は `involution <list>` で、その成分の組の一覧の識別子を付ける（例えば `involution 16-16d8`。`certificates/README.md` の Identifiers の節を見よ）。型 3-(z, f)-s の位数 3 の元は `odd prime 3-<z>-<f>-<s>`。奇素数位数の元のうち、データベースが型を記録していないものは `odd prime (type not recorded)`（S_gen の 174 個の格子）。位数 2^a のエルミートの項目は `order 2^a with g^(2^(a-1)) = -1`（a = 2, 3, 4, 5）。S_CQ については位数 4 の項目だけを挙げる。そのような元を持つ S_CQ の 238 個の格子の、位数 8、16、32 の項目は記録していない |
| | \|Aut L\| > 10^7 である 10 個の格子では、4 つの `group` フィールドはすべて空である（`null`、`false`）: これらの群はこうした同定には大きすぎ、論文の §6.3 と第 2 の論文で記述する |
| `frame` | S_(C_2^2) について、対合 σ のフレームシェイプ: その特性多項式を ∏ (x^k − 1)^(a_k) の形に書いたもの。例えば (x² − 1)^24 (x − 1)^−16 は `2^24 1^-16`、(x² − 1)^16 は `2^16`（指数が 0 の因子は書かない） |

Gram 行列と自己同型群の生成元は、リリースの添付ファイル `database-v2.tar.gz` にある。その
`index.jsonl.gz` は「ラベルの順」に並ぶが、これは数としての順である。すなわち、まず \|Aut L\| の順、次に i の順で、どちらも数として比べる（したがって 32.8.2
は 32.8.10 より前に来て、32.16.1 はすべてのラベル 32.8.i より後に来る）。

## largest-groups.jsonl

`largest-groups.jsonl` は 1 行に 1 個の JSON オブジェクトを持ち、第 2 の論文の 99 個の格子のそれぞれに 1 行が対応する。
これらは、S_gen の格子のうち |Aut L| > 2000 または |Aut L| = 512, 1024, 1536 であるものである
（それ以外の格子の群は Small Groups Library の `IdGroup` で同定される）。行は表 A1 の順に並ぶ: |Aut L| の順、次にラベルの添字 i の順である。
`verify/verify_p2.py` がこのファイルを検査する。その Gram 行列はデータベースのものである。生成元は表示の語が参照するものであり、6 個の格子（32.4608.2, 32.5760.1, 32.12288.1, 32.23328.1, 32.196608.2, 32.839808.1）では、二つの一覧が異なる。そのうち 2 個ではデータベースのほうが行列 1 個だけ多く、残りの 4 個では、それぞれの一覧に、他方に無い行列が 1 個ずつある。なお、どちらの一覧も Aut L を生成する。

| フィールド | 意味 |
|---|---|
| `label` | `32.<\|Aut L\|>.<i>`。データベース [D] のラベルであり、`classified-lattices.jsonl` のラベルと同じである |
| `id` | L の構成の識別子: `certificates/README.md`（Identifiers の節）にあるとおりの一覧の項目 `<list>:<index>`、または `known:<family>:<name>`。後者は、既知の格子（第 1 の論文の付録 A）を、出版された記述（符号、またはエルミート構造もしくは四元数的な構造）から構成したものの名前である |
| `aut` | Aut L の位数（整数） |
| `gram` | L の Gram 行列 G（32 × 32 の整数） |
| `generators_convention` | `columns: x -> g x, g^T G g = G` |
| `generators` | Aut L の生成元 g。列の規約 g^T G g = G による（第 1 の論文の §7.1 のもの） |
| `table_a1` | 表 A1 の列のうち、ラベルと位数を除くもの: `aut_factor`、`fstar`、`g_mod_fstar`、`chief`、`presentation`、`origin`。`paper/p2-lattices-99.json` と同じである |
| `presentation` | Aut L の表示。その位数は §2.3 の補題によって証明される（下記） |
| `profile` | 第 2 の論文のプロファイル（§2.2）のうち、データであるもの: `structure`（F^*(G)、G/F^*(G)、中心、導来群の位数、アーベル不変量、主列の位数、および上から順に各主因子についてその位数と `kind`（`abelian` の場合は `p`、`dim`、作用の像の位数、絶対既約性、分裂、補群の類の数。`non-abelian` の場合は `simple` と `copies`）、分裂を計算していない場合は `split_note`）。`action_on_L`（各共役類について: 位数、大きさ、中心化群の位数、フレームシェイプ、トレース、および固定格子の階数、行列式、判別式群、最小ノルム（階数 0 のときは `null`）。さらに 146880 本の最短ベクトルの上の軌道の長さ）。`construction`（小さい類の対合 σ について: 類の大きさ、C(σ) の位数、σ ~ −σ であるかどうか、F = L^σ と K = L^(−σ) の不変量と Gram 行列、それらの群の位数、C(σ) の像。最大の 2 個の群については `null` である（第 2 の論文の §4.5 がこれらを除いている。\|Aut L\| > 2 · 10^10）。`known:zeta17:s1-1` は対合が −1 だけであり、`note` がそのことを記す） |

**表示.** 表示の行列は行の規約による: 行列 S は行ベクトルに v ↦ v S で作用し、
S G S^T = G を満たす。S は、列の規約での自己同型の転置である。この規約では、
群の元の積は行列の積であり、語は因子を左から右へ掛けた積として評価する。
生成元には 0 から番号を付ける。ただし GAP の文字表現（下記）は例外である。すべての表示に共通のフィールドは
次のとおりである:

| フィールド | 意味 |
|---|---|
| `kind` | `pc`（補題 2.1）、`radical`（補題 2.3。Q = G/R に対するブロックの段が 1 つあり、R の巡回の段が続く。R は可解根基）、または `chain`（補題 2.3。段は主列に沿う） |
| `lemma` | `2.1` または `2.3` |
| `convention` | 行の規約を言葉で記したもの |
| `generator_words` | このファイルの `generators` の各行列 g（6 個の格子ではデータベースのものと同じ組ではない。上を見よ）について、同じ順に、表示の生成元の語であって積が g^T になるもの |

`kind` = `pc`:

| フィールド | 意味 |
|---|---|
| `generators` | 行列 x_0, …, x_(m−1)。主列に沿った pc 列である |
| `relative_orders` | r_0, …, r_(m−1)。その積は \|Aut L\| である |
| `relations` | x_i^(r_i) = word に対して `{"kind": "power", "i", "r", "word"}`、x_i^(−1) x_j x_i = word（i < j）に対して `{"kind": "conjugate", "i", "j", "word"}`。語は指数ベクトル (e_0, …, e_(m−1)) であり、x_0^(e_0) x_1^(e_1) ⋯ x_(m−1)^(e_(m−1)) を表す。k ≤ i では e_k = 0 である |
| `generator_words` | `relations` と同じ形の指数ベクトル |

`kind` = `radical`:

| フィールド | 意味 |
|---|---|
| `radical_generators`, `radical_relative_orders`, `radical_relations` | R の pc 表示。フィールドは `kind` = `pc` のものと同じである。相対位数の積は \|R\| である |
| `lifts` | 行列 y_1, …, y_n: Q の表示の生成元の G への持ち上げ |
| `lift_conjugates` | y^(−1) x_j y = word に対して `{"lift": a, "j", "word"}`。y は添字 a（0 から番号を付ける）の持ち上げであり、語は x についての指数ベクトルである |
| `quotient_relators` | Q の表示の関係子を GAP の文字表現で記したもの: 文字 a > 0 は y_a、−a は y_a^(−1) を表す（1 から番号を付ける） |
| `quotient_relator_values` | 各関係子について、持ち上げの上でのその値を x についての指数ベクトルで記したもの |
| `quotient_order` | \|Q\| = \|Aut L\| / \|R\| |
| `quotient_presentation_order` | `quotient_relators` が定める群の位数。`method` を伴う: 自明な部分群の上の剰余類の数え上げ、または `12-20-c63:236` については X_0 = {y}、R_0 = {y^r} とした補題 2.4（`lemma_2_4`: 生成元 y、関係子 y^r、および剰余類の数え上げで数えた ⟨y⟩ の指数）。この位数は入力である（§2.5） |
| `generator_words` | `{"radical": e, "quotient": w}`: g^T は、x についての指数ベクトル e と y についての文字の語 w の、この順の積である |

`kind` = `chain`:

| フィールド | 意味 |
|---|---|
| `generators` | すべての行列。段ごとに順に並べる |
| `stages` | `{"kind": "cyclic", "r", "generators": [a]}`、または `{"kind": "block", "k", "T_order", "components"}`。後者は、非アーベルな主因子 T^k の k 個の極小正規部分群のそれぞれに成分を 1 つずつ持つ。成分は、その `generators`（添字）、その `relators`（文字表現。成分の中で 1 から番号を付ける）、および `presentation_order` を持つ: これはその関係子が定める群の位数であり、`method` を伴う。方法は、自明な部分群の上の剰余類の数え上げか、`lemma_2_4_chain` に記した補題 2.4 の鎖である（各段階で、生成元の部分集合、その文字で書ける関係子の数、次の部分集合、剰余類の数え上げで数えた指数。最後の段階では、剰余類の数え上げで数えた位数）。これらの位数は入力である（§2.5） |
| `relations` | x_a^r = word に対して `{"kind": "power", "stage", "a", "r", "word"}`。成分の関係子をその生成元の上で評価したもの = word に対して `{"kind": "relator", "stage", "component", "relator", "word"}`。x_a^(−1) x_b^(−1) x_a x_b = word（a と b は異なる成分に属する）に対して `{"kind": "commutator", "stage", "a", "b", "word"}`。x_a^(−1) x_b x_a = word（b は後の段に属する）に対して `{"kind": "conjugate", "stage", "a", "b", "word"}`。語は組 `[generator, exponent]` のリストであり、冪を左から右へ掛けた積を表す。語に現れる生成元はすべて `stage` より後の段に属する |
| `order_bound` | 補題 2.3 の上界: 巡回の段についての r の積と、ブロックの段についての \|T\|^k の積。これは \|Aut L\| に等しい |
| `generator_words` | `relations` と同じ形の組のリスト |

関係の右辺と生成元の語は GAP で求めた（多重巡回の指数、
有限表示群への準同型、段を順にたどるふるい分け）。`verify/verify_p2.py` はそれらを整数行列で評価するので、
誤りがあれば検出される。求め方は結論に
影響しない。

## データの読み方

Python の標準ライブラリだけを用い、リポジトリのルートから、`largest-groups.jsonl` の格子を 1 個読み、
その Gram 行列が偶で対称であること、およびすべての生成元が等長写像であることを検査する。

```python
import json
rows = {r['label']: r for r in map(json.loads, open('data/largest-groups.jsonl'))}
r = rows['32.512.1']
G, gens = r['gram'], r['generators']            # 32 x 32 integers; generators act on columns: x -> g x
n = len(G)
mul = lambda A, B: [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
tr = lambda A: [list(c) for c in zip(*A)]
assert all(G[i][j] == G[j][i] for i in range(n) for j in range(n)) and all(G[i][i] % 2 == 0 for i in range(n))
assert all(mul(mul(tr(g), G), g) == G for g in gens)   # each generator is an isometry of L
print(r['aut'], len(gens))
```

## リリースの添付ファイルのフィールド

添付ファイル `database-v2.tar.gz`、`general-part-records.tar.gz`、`general-part-groups.tar` は、それぞれの README を含む。
それ以外については次のとおりである。

`type-8-24-classes.tar` は 4 個のファイル `classes-0.jsonl.gz`, ..., `classes-3.jsonl.gz` からなり、P = A_1^24 についての
良い部分空間 W の 376306 個の軌道を収める（第 1 の論文の付録 B.3）。
軌道 i はファイル `classes-<i mod 4>` にある。

| フィールド | 内容 |
|---|---|
| `i` | 軌道の番号 0, ..., 376305 |
| `stab` | Aut P の中での \|Stab(W)\| |
| `gram` | K_W の Gram 行列（階数 24、偶、行列式 256、最小ノルム 4）。LLL 簡約してある |
| `aut` | \|Aut K_W\| |
| `index` | [N : N']。ここで N = √2 K_W^# であり、N' は K_W^# のノルム 2 の 720 本のベクトルで張られる |
| `method`, `cert` | K_W を同定した方法: `graph`（720 本のベクトルのグラフの標準形。`cert` はそのハッシュ）、`orbit`（同じものに N / N' の軌道を加えたもの。`cert` はその最小の元のハッシュを加える）、または `qfauto`（`cert` は不変量のハッシュであり、等長性は `qfisom` で判定した） |
| `orbit` | そのグラフの自己同型群の中での Aut K_W の指数（`qfauto` の場合は計算しておらず、0） |
| `census` | K_W^# / K_W の 270 個の極大完全特異部分空間 T' について、K_W + T' のルートの本数ごとの T' の個数 |
| `np` | K_W^# / K_W の 135 個の 0 でない特異な類 p について、p の中のノルム 2 のベクトルの本数ごとの類の個数 |

`gaussian-classes.jsonl.gz` は、命題 4.6 のルートを持たないエルミート Z[i]-格子 76742 個のそれぞれに 1 行を持つ。

| フィールド | 内容 |
|---|---|
| `key` | 第 1 の論文の §7.2 のキー: `k` と 16 進 24 桁（離散的な組の標準形の Gram 行列）、または `D` と 16 進 16 桁（不変量 D） |
| `D` | 不変量 D（16 進 16 桁） |
| `gram`, `W` | トレース格子の Gram 行列 G と、i 倍の行列 W。行の規約による: W G W^T = G、W^2 = −1 |
| `U` | ユニタリ群の位数 |
| `antilinear` | その格子が複素共役と等長であるかどうか |
| `f`, `w` | `antilinear` ならば f = 1、そうでなければ 2。w は命題 4.6 の重み w = f / \|U\| |

`gaussian-walk-grams-1.jsonl.gz` と `gaussian-walk-grams-2.jsonl.gz` は、2 本の近傍のランダムウォークが出会った類について、
上と同じフィールド `key`、`gram`、`W` を持ち、さらに `id` を持つ: 第 1 のウォークでは `k11-<n>`、第 2 のウォークでは `k12-<n>`
である（11 と 12 はウォークの種）。ここで n は、その類に出会った時点でウォークが訪れていたルートを持たない格子の数である。
