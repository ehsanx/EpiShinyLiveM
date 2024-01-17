# Define the paths
getwd()
app_path <- paste0(getwd(), "/app/")
output_path <- paste0(getwd(), "/docs")

# Use the paths in the export function
shinylive::export(app_path, output_path)
httpuv::runStaticServer(paste0(getwd(), "/docs"))
servr::httd(paste0(getwd(), "/docs"))

