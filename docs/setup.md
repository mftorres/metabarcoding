# Ubuntu/WSL setup and obtaining the project

Use Ubuntu for the commands below. If Ubuntu is not installed, follow
[Microsoft's WSL installation instructions](https://learn.microsoft.com/windows/wsl/install).

## Access Windows folders from Ubuntu

WSL exposes Windows drive `C:\` at `/mnt/c/`. Replace `USERNAME` with your
**Windows account folder name**, which may differ from your Ubuntu username.
For example, Windows `C:\Users\USERNAME\Downloads` is available in Ubuntu as:

```bash
cd /mnt/c/Users/USERNAME/Downloads
```

Your Windows Documents folder is normally `/mnt/c/Users/USERNAME/Documents`.
If Windows redirects Documents into OneDrive, use the actual location, such as
`/mnt/c/Users/USERNAME/OneDrive/Documents`. Use forward slashes in Ubuntu and
quote the entire path when a folder name contains spaces. Run `ls /mnt/c/Users`
to find the Windows account folder, and `explorer.exe .` to open the current
Ubuntu location in Windows Explorer. Ubuntu's `~` refers to your Linux home,
not your Windows user folder.

## Download and install Linux Miniconda

Use the [official Miniconda installer directory](https://repo.anaconda.com/miniconda/).
For this project's Linux x86_64 / WSL2 environment, download
[Miniconda3-latest-Linux-x86_64.sh](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)
with your browser and save it in Windows Downloads. Choose the **Linux `.sh`**
installer, even though your computer runs Windows. Check `uname -m` in Ubuntu;
these commands and the environment recipe target `x86_64`.

If `conda --version` already works in Ubuntu, use that installation instead of
installing another Conda. Otherwise, open Ubuntu and run:

```bash
sudo apt update && sudo apt install -y curl ca-certificates bzip2 util-linux git
cd /mnt/c/Users/USERNAME/Downloads
```

Alternatively, download directly from Ubuntu in that Windows Downloads folder:

```bash
curl -fLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sha256sum Miniconda3-latest-Linux-x86_64.sh
```

Compare the checksum with the SHA256 shown beside the same filename in the
official installer directory. Continue only when it matches. Run the installer:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow its prompts, install under your Linux home (normally `~/miniconda3`),
and allow Bash initialisation. Close and reopen Ubuntu, then run
`conda --version`. If necessary, run `~/miniconda3/bin/conda init bash` and reopen
Ubuntu. Installing Conda does not download the analysis scripts or databases.

## Clone the metabarcoding repository

The project is at [github.com/mftorres/metabarcoding](https://github.com/mftorres/metabarcoding).
Install Git with `sudo apt install -y git` if needed. Clone into your Windows
Documents folder; `devfull` contains the two-module workflow documented here:

```bash
cd /mnt/c/Users/USERNAME/Documents
git clone --branch devfull https://github.com/mftorres/metabarcoding.git
cd metabarcoding
pwd
ls code config data
```

Use your actual OneDrive Documents path if applicable. Run `git clone` once, into
a location where `metabarcoding` does not already exist. An existing clone can be
opened directly with `cd /mnt/c/Users/USERNAME/Documents/metabarcoding`.
A browser download via GitHub's **Code → Download ZIP** is another option; select
`devfull` first and extract the ZIP before using the scripts.

## Performance: run analyses on the Linux filesystem

Under WSL, running analyses in the local Linux filesystem, such as
`~/metabarcoding`, may be faster than working in Windows-mounted folders such as
`/mnt/c/Users/USERNAME/Documents/metabarcoding`, especially when reading and
writing many files. This follows [Microsoft's WSL filesystem guidance](https://learn.microsoft.com/windows/wsl/filesystems).
For a new project in your Linux home, use:

```bash
cd ~
git clone --branch devfull https://github.com/mftorres/metabarcoding.git
cd ~/metabarcoding
```

Place the input reads and reference databases in this Linux project too, and
write analysis outputs there to benefit from local Linux I/O. You can still
access the folder from Windows by running `explorer.exe .` in Ubuntu. Performance
varies with the workload and storage hardware. Choose the location before
starting analyses: existing manifests contain absolute paths, so moving a
completed analysis may require recreating its post-classification outputs.

From the project root, create and activate the environment:

```bash
conda config --set channel_priority strict
conda env create -f config/taxonomy-environment.yml
conda activate taxonomy
```

Raw sequence data and large reference databases are supplied separately from Git.
Follow the [MIDORI instructions](../data/sintax/README.md) for Read processing and the
[BOLD instructions](../data/bold/README.md) for Postprocessing's `--bold` mode.
The [offline HTML guide](metabarcoding-guide.html) includes the full setup and
workflow instructions.
