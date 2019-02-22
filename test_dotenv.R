require(dotenv)
dotenv::load_dot_env(file=".env") #change the file name from the default if you want to use a different file in a script

#load the environment parameters that you care about into the relevant variables
igraph_path <- Sys.getenv("IGRAPH_PATH")
r_architecture <- Sys.getenv("R_ARCH")

#here's all the variables in the environment
Sys.getenv()
