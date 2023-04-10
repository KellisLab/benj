-- module use /path/to/benj/modules; module load miniconda
-- shamelessly taken from https://github.com/AaltoSciComp/lmod-conda/blob/master/conda.lua

-- The base path to the conda environment
local root = "~/data/miniconda3"

-- only one conda module can be loaded at once:
family("conda")



--

if myShellName() == "csh" then
   execute{cmd=table.concat({"source ", pathJoin(root, "bin/activate.csh")}),
           modeA={"load"}}

else
   -- The extra '' is needed, or else $@ is set to "module load" and then
   -- "conda activate" gives an error message.
   execute{cmd=table.concat({"source ", pathJoin(root, "bin/activate ''")}),
           modeA={"load"}}

end

execute{cmd="conda deactivate", modeA={"unload"}}
execute{cmd="unset conda", modeA={"unload"}}
execute{cmd=table.concat({"echo ",  pathJoin(root, "condabin")}), modeA={"unload"}}


-- TODO: remove condabin from PATH, the following does NOT work
-- because "conda deactivate" runs later and this is what adds it.
remove_path("PATH", pathJoin(root, "condabin"))
