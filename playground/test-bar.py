'''
author:        Wang Chenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-10-02
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''
import sys
import time
from multiprocessing import Pool


def draw_bar(done, total, elapsed, width=40, prefix="进度"):
    """在终端同一行刷新进度条（纯标准库）"""
    ratio = done / total if total else 1.0
    filled = int(width * ratio)
    bar = "█" * filled + "-" * (width - filled)

    speed = done / elapsed if elapsed > 0 else 0.0          # 每秒完成数
    eta = (total - done) / speed if speed > 0 else 0.0      # 预计剩余秒数

    line = (f"\r{prefix}: |{bar}| {done}/{total} "
            f"({ratio:5.1%}) {speed:6.1f}it/s ETA {eta:5.1f}s")
    # 末尾补空格，防止残留上一次更长的内容
    sys.stderr.write(line.ljust(100))
    sys.stderr.flush()


def work(x):
    """模拟一个耗时任务"""
    time.sleep(0.05)
    return x * x


def main():
    total = 100
    tasks = list(range(total))

    start = time.time()
    done = 0
    results = []

    with Pool(8) as pool:
        # imap_unordered：谁先完成谁先返回，进度条才平滑
        for r in pool.imap_unordered(work, tasks, chunksize=1):
            results.append(r)
            done += 1
            draw_bar(done, total, time.time() - start)

    sys.stderr.write("\n")
    print("结果数量:", len(results), "前5个:", results[:5])


if __name__ == "__main__":
    main()