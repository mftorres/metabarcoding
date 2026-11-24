# Setup: Ubuntu (WSL), Miniconda and the project

You need to do this once per computer. All commands are typed in the **Ubuntu**
window (WSL, the Windows Subsystem for Linux). If Ubuntu is not installed, follow
[Microsoft's WSL installation instructions](https://learn.microsoft.com/windows/wsl/install).

## 1. Find your Windows folders from Ubuntu

Ubuntu sees the Windows `C:` drive as `/mnt/c/`. Replace `USERNAME` below with your
**Windows account folder name** (run `ls /mnt/c/Users` to see it). For example,
Windows `C:\Users\USERNAME\Downloads` is:

```bash
cd /mnt/c/Users/USERNAME/Downloads
```

- If Windows keeps Documents in OneDrive, the path is `/mnt/c/Users/USERNAME/OneDrive/Documents`.
- Use forward slashes `/`, and put quotes around a path that contains spaces.
- `~` means your **Linux** home folder, not your Windows user folder.
- `explorer.exe .` opens the current Ubuntu folder in Windows Explorer.

## 2. Install Miniconda (the software manager)

If `conda --version` already prints a version, skip this step.

```bash
sudo apt update && sudo apt install -y curl ca-certificates bzip2 util-linux git
cd ~
curl -fLO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sha256sum Miniconda3-latest-Linux-x86_64.sh
```

Compare the printed checksum with the one listed beside that file on the
[Miniconda download page](https://repo.anaconda.com/miniconda/); continue only if they
match. Then run the installer, accept the default location (`~/miniconda3`) and say
yes to initialising Bash:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Close and reopen Ubuntu, then check `conda --version`.

## 3. Get the project

Choose where the project lives **before** you start analysing; moving it later is
possible but large data folders are slow to move.

Analyses run faster inside the Linux file system (`~/…`) than in Windows folders
(`/mnt/c/…`), because WSL reads Windows files more slowly
([Microsoft's explanation](https://learn.microsoft.com/windows/wsl/filesystems)).
Recommended:

```bash
cd ~
git clone --branch devfull https://github.com/mftorres/metabarcoding.git
cd ~/metabarcoding
ls code config data
```

You should see the `code`, `config` and `data` folders. Every command in these guides
is run from this project folder. To open it in Windows Explorer, run `explorer.exe .`.
(To keep it in Windows Documents instead, `cd /mnt/c/Users/USERNAME/Documents` before
`git clone`.)

## 4. Create the two software environments

An environment is a self-contained set of programs. The project uses two:

| Environment | Contains | Used for |
| --- | --- | --- |
| `taxonomy` | fastplong, VSEARCH, SeqKit, MAFFT, kc-align, ClipKIT, IQ-TREE, Python libraries | everything you run yourself |
| `cutadapt` | cutadapt | primer finding and trimming in Module 1; the pipeline calls it automatically |

cutadapt has its own environment because it is not yet available for the Python
version used by `taxonomy`. You never need to activate `cutadapt` yourself.

```bash
cd ~/metabarcoding
conda config --set channel_priority strict
conda env create -f config/taxonomy-environment.yml
conda env create -f config/cutadapt-environment.yml
conda activate taxonomy
fastplong --version && vsearch --version 2>&1 | head -1 && conda run -n cutadapt cutadapt --version
```

The last line should print three version numbers. If an environment already exists,
`conda env create` stops with an error; that is fine, use the existing one.

## 5. Download the databases

Large databases are not stored in Git. Follow:

- [MIDORI (needed by Module 1)](../data/sintax/README.md)
- [BOLD (needed by Module 2 with `--bold`)](../data/bold/README.md)

Then continue with [Module 1 · Read processing](read_processing.md).
