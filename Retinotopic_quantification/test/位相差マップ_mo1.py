import os
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt

#-----------------------------------
# 1. TIFFファイルを読み込み、チャネルを抽出
#-----------------------------------
tiff_path = "./test4.tif"
if not os.path.exists(tiff_path):
    raise FileNotFoundError(f"TIFF file not found: {tiff_path}")

image = tiff.imread(tiff_path)

if len(image.shape) < 2:
    raise ValueError("Invalid TIFF file: Image must have at least two dimensions.")

num_channels = image.shape[0] if len(image.shape) == 3 else 1
if num_channels < 3:
    raise ValueError("The TIFF file does not contain at least 3 channels.")

# チャネル1（最初のチャネル）とチャネル3（3番目のチャネル）を抽出
channel1_orig = image[0, ...].astype(float)
channel3_orig = image[2, ...].astype(float)

# 各チャネルを[0,1]の範囲に正規化
channel1_orig = (channel1_orig - channel1_orig.min()) / (channel1_orig.max() - channel1_orig.min())
channel3_orig = (channel3_orig - channel3_orig.min()) / (channel3_orig.max() - channel3_orig.min())

#-----------------------------------
# 2. エッジアーティファクトを除去するために画像をパディング
#    2の累乗サイズのキャンバスを使用し、平均値でパディング
#-----------------------------------
# 元の画像サイズを取得
orig_height, orig_width = channel1_orig.shape

# キャンバスサイズを決定：最大次元の1.3倍以上で最も近い2の累乗
minFFTsize = 1.3 * max(orig_width, orig_height)
canvas_size = 1
while canvas_size < minFFTsize:
    canvas_size *= 2

# 平均値でパディングする関数

def pad_with_mean(image, canvas_size):
    # image: 2D配列, canvas_size: int（正方形キャンバスを想定）
    orig_h, orig_w = image.shape
    padded = np.empty((canvas_size, canvas_size), dtype=image.dtype)
    
    # 元の画像の平均値を計算し、埋める
    image_mean = np.mean(image)
    padded.fill(image_mean)
    
    # 元の画像を中央に配置
    start_y = (canvas_size - orig_h) // 2
    start_x = (canvas_size - orig_w) // 2
    padded[start_y:start_y+orig_h, start_x:start_x+orig_w] = image
    return padded

# 両チャネルをパディング
channel1 = pad_with_mean(channel1_orig, canvas_size)
channel3 = pad_with_mean(channel3_orig, canvas_size)

print(f"Original dimensions: {orig_width}x{orig_height}, Padded canvas size: {canvas_size}x{canvas_size}")

#-----------------------------------
# 3. パディングした画像にフーリエ変換を適用
#-----------------------------------
F1 = np.fft.fft2(channel1)
F3 = np.fft.fft2(channel3)


#-----------------------------------
# 3.5. 100-200µmの波長に焦点を当てる周波数マスクを作成
#-----------------------------------
# 与えられた値: 1 µm = 0.6402 ピクセル (元の画像)
# したがって：
#   100 µm ≈ 100 * 0.6402 ≈ 64 ピクセル
#   200 µm ≈ 200 * 0.6402 ≈ 128 ピクセル

f_lower = 1 / 128  # 周波数（ピクセル単位の周期）
f_upper = 1 / 64   # 周波数（ピクセル単位の周期）

# パディングされた画像の周波数グリッドを作成
freqs = np.fft.fftfreq(canvas_size, d=1)  # d=1（ピクセル間隔）
freqs = np.fft.fftshift(freqs)            # ゼロ周波数を中心にシフト
fx, fy = np.meshgrid(freqs, freqs)
radial_freq = np.sqrt(fx**2 + fy**2)

# バンドパス マスクを作成：周波数が f_lower から f_upper の範囲にある場合に True
band_mask = (radial_freq >= f_lower) & (radial_freq <= f_upper)

# フーリエ変換にマスクを適用
F1_band = F1 * band_mask
F3_band = F3 * band_mask

#-----------------------------------
# 4. 位相解析（バンドパスフィルタ適用後）
#-----------------------------------
# バンド制限データの位相を計算
phase1_band = np.angle(F1_band)
phase3_band = np.angle(F3_band)

# 選択したバンドの位相差を計算（[-π, π] にラップ）
phase_diff_band = np.angle(np.exp(1j * (phase3_band - phase1_band)))
phase_diff_band_shifted = np.fft.fftshift(phase_diff_band)

# バンド制限データのクロスパワースペクトルを計算
epsilon = 1e-8  # ゼロ除算を避けるための小さい値
cross_power_band = F1_band * np.conj(F3_band) #複素数の共役複素数を計算する
cross_power_band_norm = cross_power_band / (np.abs(cross_power_band) + epsilon) #正規化し位相情報のみを残す

# バンドの位相相関マップを計算
corr_map_band = np.fft.ifft2(cross_power_band_norm)
corr_map_band_abs = np.abs(corr_map_band)
corr_map_band_shifted = np.fft.fftshift(corr_map_band_abs)

#-----------------------------------
# ５. プロット描画
#-----------------------------------
plt.figure(figsize=(14, 4))

plt.subplot(1, 3, 3)
plt.title("Band-Pass Phase Difference (Ch3 - Ch1)")
plt.xlabel("X (pixels)")
plt.ylabel("Y (pixels)")
plt.imshow(phase_diff_band_shifted, cmap='twilight', origin='lower')
plt.colorbar(label="Phase Difference (radians)")

plt.tight_layout()
plt.show()
