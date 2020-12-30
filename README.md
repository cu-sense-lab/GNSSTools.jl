# GNSSTools

A general Julia package for working with GNSS data. Support course acquisition for L5Q signals. New will signals will be added on an as needed basis. Local signal generation and course acquisition methods support multithreading. See instructions for installing this package.

**NOTE:** This is a private repo, but you should still have access to it if logged in. Contact me if you don't.

# Installation

## Allowing Julia to Access SeNSe GitLab Repos

Since this is a private repo, the installation instructions are different than for publically registered Julia packages. Add the following to your `.bashrc` file.

```bash
# Startup ssh-agent and add RSA key for Julia's Package Manager
eval `ssh-agent -s`
ssh-add ~/.ssh/id_rsa  # This could also be the path to your RSA key if it is named differently
```

The above ensures that Julia's package manager can access SeNSe's internal GitLab using the same ssh key you use to make commits. Once you've done this, restart your terminal session by logging out and back in. You will now be able to install internal packages directly through Julia's package manager.

**NOTE:** If your ssh key for GitLab changes, make sure to update the above to point to it.

## Installing GNSSTools

From terminal type `julia` and the type `]` to enter the package manager. Then copy/paste the following:

```julia
add git@192.168.3.66:bilardis/GNSSTools.jl.git
```

Julia will automatically install the neccessary dependencies and GNSSTools. Once complete, you can return the Julia REPL by hitting the `backspace` key. To update GNSSTools if newer versions are released, enter the package manager again by typing `]` and use `update GNSSTools`.
