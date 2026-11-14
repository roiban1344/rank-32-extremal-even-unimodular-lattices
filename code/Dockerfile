# The environment in which the programs were run, as a container (linux-64).
# It installs exactly the packages of environment.lock.txt (conda-forge, with MD5 digests) on a pinned
# micromamba image; see README.md, "Requirements".
#
# From the repository root (on machines other than linux-64, add --platform linux/amd64 to both):
#
#   docker build -t rank32 code
#   docker run --rm -v "$PWD":/work:ro -w /work rank32 python -m pytest -q -p no:cacheprovider code/tests
#
# プログラムを走らせた環境を、コンテナにしたもの（linux-64）。版を固定した micromamba のイメージの上に、
# environment.lock.txt のパッケージ（conda-forge、MD5 ダイジェスト付き）をそのまま入れる。README.ja.md の
# 「必要なもの」を参照。コマンドはリポジトリの一番上で走らせる（linux-64 以外のマシンでは、両方に
# --platform linux/amd64 を付ける）。
FROM mambaorg/micromamba:2.10.0-ubuntu24.04@sha256:c3897a06f87c1c13f038188b3665a8c4162da62eac5109ae6a939d1e290dcefa

# The lock file stays in the image, at /opt/rank32/environment.lock.txt, as the list of what it holds.
# ロックファイルは、イメージの中身の一覧として /opt/rank32/environment.lock.txt に残す。
COPY environment.lock.txt /opt/rank32/environment.lock.txt
RUN micromamba install -y -n base -f /opt/rank32/environment.lock.txt \
    && micromamba clean -a -y

WORKDIR /work
