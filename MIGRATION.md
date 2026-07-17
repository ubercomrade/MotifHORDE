# Migration From Python Outputs

The Julia cutover deliberately does not include a pickle/joblib compatibility
loader. Python model objects are runtime-specific and are not a stable
interchange format.

To use an old result, regenerate the model from the original discovery output,
export a MEME PFM, or rerun the discovery tool with the Julia CLI. New model
bundles are directories produced by `Mimosa.writemodel`; they are not renamed
`.pkl` files.

SiteGA training is now owned by an independent Julia project. MotifHORDE only
invokes its versioned process contract and reads the listed Mimosa-compatible
models from its manifest.
