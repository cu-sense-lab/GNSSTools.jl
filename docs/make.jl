using Documenter
using DocumenterMarkdown
using GNSSTools


"""
    fix_eq_format(file_name)

Fix the format of the equations in the Markdown
files to be compatible with the GitLab equation
format requirement.
"""
function fix_eq_format(file_name)
    lines = readlines(file_name)
    for i in 1:length(lines)
        if occursin("\$\$", lines[i])
            lines[i] = "```math\n"
            i += 1
            while ~occursin("\$\$", lines[i])
                i += 1
            end
            lines[i] = "\n```"
        end
    end
    file = open(file_name, "w")
    for line in lines
        if line == ""
            write(file, "\n")
        else
            write(file, string(line, "\n"))
        end
    end
    close(file)
end


"""
    fix_format_in_files(directory, folders)

Fix the format of equations in each Markdown in
folders within `directory`.
"""
function fix_format_in_files(directory, folders)
    if isa(folders, Tuple)
        folers = [folders]
    end
    for folder in folders
        src = string(directory, "/", folder[1])
        dest = string(directory, "/", folder[2])
        files = readdir(src)
        for file in files
            if endswith(file, ".md")
                fix_eq_format(string(dest, "/", file))
            end
        end
    end
end

build_name = "GNSSTools.jl.wiki"

# Generate API and other page Markdown files
makedocs(sitename="GNSSTools Documentation",
         # format=Documenter.HTML(prettyurls=true)
         format=Markdown(),
         build=build_name,
         clean=false)

# Correct equation format in Markdown files
@info "Correcting equations."
directory = string(pwd(), "/docs")
folders = [("src", build_name),
           ("src/API", string(build_name, "/API"))]
fix_format_in_files(directory, folders)
