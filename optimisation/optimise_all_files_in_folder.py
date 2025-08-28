import os
import subprocess
import cv2
import numpy as np


def execute_script_for_png_files(folder_path, config_path, results_folder_base):
    png_files = [f for f in os.listdir(folder_path) if f.endswith('.png')]

    for png_file in png_files:
        name = os.path.splitext(png_file)[0]

        img = cv2.imread(f'{folder_path}/{name}.png', cv2.IMREAD_UNCHANGED)
        change = False
        image_size = img.shape[0]

        if np.iinfo(img.dtype).max == 255:
            print(f'{name} 8 bits image will not be processed.')
            continue
        else:
            print(f'{name} will be processed!')

        if img.shape[0] != img.shape[1]:
            m = min(img.shape[0], img.shape[1])
            img = img[:m, :m]
            change = True
            image_size = m
        if img.shape[0] > 1024:
            img = img[:1024, :1024]
            change = True
            image_size = 1024
        if len(img.shape) == 3:
            img = img[:, :, 0]
            change = True

        if change:
            name = f'{name}_crop'
            cv2.imwrite(f'{folder_path}/used/{name}.png', img)

        command = [
            "python", "optimise.py",
            f"--config_path={config_path}",
            f"--optim.image_file_name={folder_path}/used/{name}.png",
            f"--log.results_folder={results_folder_base}{name}",
            # f"--optim.image_size={image_size}",
            f"--multiscale.samples_per_res=[8000, 10000, 10000, 30000, 50000]"
        ]

        print(f"Executing command for {png_file}: {' '.join(command)}")
        subprocess.run(command)


if __name__ == "__main__":
    # PNG files
    folder_path = "data/dem/SRTM"

    config_path = "configs/default.yml"
    results_folder_base = "results/gsplat/dem/SRTM/timing/start108k/"

    execute_script_for_png_files(folder_path, config_path, results_folder_base)
