# 階数 32 の極値偶ユニモジュラー格子

[English](README.md) | 日本語
<!-- sync-extra: 32 32 （二つの論文の表題の日本語訳にある数） -->

このリポジトリは二つの論文に付属する。二つの論文は、どの版でも一緒に公開する。

- *The extremal even unimodular lattices of rank 32: their number, and the classification of those with
  automorphisms other than ±1*（階数 32 の極値偶ユニモジュラー格子: その個数と、±1 以外の自己同型を持つものの分類。`paper/extremal-lattices-rank-32.pdf`、第 2 の論文では [P1]）
- *Automorphism groups of extremal even unimodular lattices of rank 32*（階数 32 の極値偶ユニモジュラー格子の自己同型群。`paper/automorphism-groups-rank-32.pdf`）。群の最も大きい 99 個の格子を扱う。

論文と補助ノートには日本語版がある（`paper/extremal-lattices-rank-32.ja.pdf`、`paper/automorphism-groups-rank-32.ja.pdf`、`paper/assignment-count.ja.pdf`）。定理・式・表・計算の番号は英語版と同じである。

**結果.** 階数 32 の極値偶ユニモジュラー格子の等長類は、ちょうど 11044268 個ある（定理 A）。そのうちちょうど 143329 個が ±1 以外の自己同型を持つ。論文はこれらの格子を分類し、それぞれの自己同型群の位数を決定する（定理 B）。証明は計算機援用である。論文の表 7 に証明が依拠する主張を、表 10 にソフトウェアについて仮定することを挙げてある。

第 2 の論文は、群の位数が 2000 を超えるか 512、1024、1536 に等しい 99 個の格子の群を記述する（それ以外の格子の群は Small Groups Library の `IdGroup` で同定される）。記述するのは、一般 Fitting 部分群、主列、位数を証明した表示、共役類、最短ベクトルの上の軌道である。なお、99 個のうち 27 個は文献にあり、残りの 72 個は、我々が読んだ文献には見つからなかった。

## 内容

| パス | 内容 |
|---|---|
| `paper/` | 二つの論文（Typst のソースと PDF）と、その計算の一覧表（`computations.typ` は第 1 の論文の表 11、`computations-p2.typ` は第 2 の論文の表 11）、第 2 の論文の表 A1 のデータ（`p2-lattices-99.json`）、および数え上げの以前の整理についての補助ノート（成分の一覧を共有する。`assignment-count.pdf`、第 1 の論文では [Aux]） |
| `certificates/` | 検査が読む certificate（検証用データ）。主張または計算ごとに 1 ファイル。`certificates/README.md` を参照 |
| `data/classified-lattices.jsonl` | ±1 以外の自己同型を持つ 143329 個の格子: ラベル、自己同型群の位数と（同定されている場合は）構造、およびその格子が属する定理 B の部分。`data/README.md` を参照 |
| `data/largest-groups.jsonl` | 第 2 の論文の 99 個の格子: ラベル、Gram 行列、自己同型群の生成元、表示、およびその論文の表 A1 の列。`data/README.md` を参照 |
| `verify/` | `verify.py` は certificate とデータを第 1 の論文と照らして検査し、`verify_p2.py` は第 2 の論文のデータを整数で検査する。`verify/README.md` を参照 |
| `LICENSE` | `verify/` のコードは MIT、文章とデータは CC BY 4.0 |
| `MANIFEST.json` | この版のすべてのファイルの SHA-256 ダイジェスト |

リリースの添付ファイル（論文の表 12）には、大きなデータと研究記録がある。

| 添付ファイル | 内容 |
|---|---|
| `research-record.tar.gz` | 一覧と certificate を作ったプログラム（変更なし）。計算ごとの索引と検証のノートが付く（論文では [Rec]） |
| `research-runs.tar.gz` | それらのプログラムの実行の記録 |
| `gaussian-classes.jsonl.gz` | 命題 4.6 の、ルートを持たないエルミート Z[i]-格子 76742 個 |
| `gaussian-walk-grams-1.jsonl.gz`, `gaussian-walk-grams-2.jsonl.gz` | Gauss 構造の一覧を作った 2 本の近傍のランダムウォークが出会った Gram 行列（フィールドは `data/README.md` にある） |
| `type-8-24-classes.tar` | Niemeier 格子 A_1^24 についての良い部分空間 W の 376306 個の軌道（第 1 の論文の付録 B.3）。それぞれに K_W の Gram 行列、\|Aut K_W\|、W の固定部分群の位数が付く |
| `general-part-groups.tar` | \|Aut L\| <= 10^7 である S_gen の 12536 個の格子の生成元、表示、共役類、群の同定。ラベルごとに収め、README が付く |
| `general-part-records.tar.gz` | S_gen の 17736 件の全自己同型群の記録（第 1 の論文の §6.3）。それぞれにラベル、\|Aut L\|、共役類、項目の行を収め、README が付く |
| `database-v2.tar.gz` | データベース [D]、版 2: S の 143329 個の格子それぞれの Gram 行列と自己同型群の生成元。README、ラベルによる索引、マニフェストが付く。`verify/verify_database.py` はこれをこのリポジトリと照らして検査し、ラベルから格子を引く |

## 検査

```
python3 verify/verify.py
python3 verify/verify_p2.py
```

どちらも、Python 3 とその標準ライブラリ（Python 3.12 と 3.14 で試した）およびこのリポジトリのファイルだけで動き、すべての検査が通ると終了ステータス 0 で終わる。`verify.py` は有限な certificate を検査し、完全性の等式をその項から足し直し、第 1 の論文の §6 の集計を保存した行から再実行し、データがその論文と一致することを検査する。なお、群・一覧・質量の計算そのものは繰り返さない（`verify/README.md` と第 1 の論文の §7.4 を参照）。
`python3 verify/verify_database.py <path to database-v2.tar.gz>` は、リリースの添付ファイルであるデータベースを検査する。そのマニフェストと索引、すべてのラベル（`data/classified-lattices.jsonl` と照合）、すべての Gram 行列（偶、正定値、行列式 1）、すべての生成元（g^T G g = G）を検査し、|Aut L| <= 10^6 については生成元が位数 |Aut L| の群を生成することも検査する。全部でおよそ 30 分かかる。また、`--label 32.<|Aut L|>.<i>` を付けると、1 個の格子の記録を表示する。
`verify_p2.py` は、与えられた生成元が 99 個の格子の自己同型であること、各表示の関係がそれらについて成り立つこと、各表示の位数が |Aut L| であることを、整数で検査する。剰余類の数え上げで得た位数は、記載の値をそのまま用いる（計算し直さない。第 2 の論文の §2.5）。

## 再現

安価なものから高価なものまで、三つの水準がある。

1. **リポジトリ.** 上のとおり、`python3 verify/verify.py` と `python3 verify/verify_p2.py`（およそ 30 秒と 1 分）。`verify/README.md` に、すべての検査と、それぞれが入力に取るものを挙げてある。
2. **一つの計算.** 論文の表 11（`paper/computations.typ`）は、各計算 C1, C2, ... について、ここにあるその certificate と、リリースの添付ファイル `research-record.tar.gz` の中の置き場 `C<n>/` を示す。そこにある `INDEX.md` は、各計算（第 1 の論文の C1, C2, ... と、第 2 の論文の A1, ..., A16。後者は `A<n>/` の下）のプログラムを挙げ、書かれている場合には再実行のノート `REPRODUCE.md` も挙げる。このノートには、コマンド、入力とそのダイジェスト、出力、再実行が再現すべきものが書いてある。なお、プログラムは実行したときのまま、研究時のパスで保存してある。実行には SageMath、PARI/GP、GAP、bliss が必要で、数時間から数日かかったものもある（表 11）。
3. **データ.** リリースの添付ファイルには、certificate を書き出した元の一覧と記録がある。また、その大きさと SHA-256 ダイジェストは `MANIFEST.json`（フィールド `release_assets`）と、第 1 の論文の表 12 にある。certificate は、フィールド `computations` に、由来する計算の名前を持つ。

## 組版

`paper/` の PDF は、隣にあるソースから Typst 0.14.2 で組んだ。たとえば `typst compile paper/extremal-lattices-rank-32.typ`。PDF は、Typst に付属するフォント New Computer Modern、New Computer Modern Math、Libertinus Serif、DejaVu Sans Mono と、別にインストールが必要な DejaVu Sans を使う。なお、Typst の別の版や別のフォントでは、文章は同じであるが、改行と改頁が変わることがある。

日本語版（`*.ja.typ`）は、これらに加えて IPAex 明朝と IPAex ゴシックを使う。これらも別にインストールが必要である。

## 版

版には `vYYYY-MM-DD.N` の形のタグを付ける。日付は、その版を作った研究用リポジトリのコミットの日付（日本標準時）であり、`MANIFEST.json` の `source_date_epoch` でもある。N は、その日のうちの版の通し番号である。後の版が前の版を訂正または拡張することがある。論文は表題の頁にその版を記す。
`MANIFEST.json` は、その版を作った研究用リポジトリのコミットも記録する。なお、研究用リポジトリは公開しておらず、コミットはビルドを識別するためだけに記す。

## 著者と利用許諾

roiban1344 が研究を指揮し、公開に責任を負う。Claude (Anthropic) と OpenAI Codex を含む AI システムが、研究、コード、検査、執筆に寄与した。

コードは MIT ライセンスで、文章とデータは CC BY 4.0 で配布する。[LICENSE](LICENSE) を参照。
