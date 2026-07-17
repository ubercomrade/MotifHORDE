"""Thin, explicit aliases around Mimosa.jl's public model API."""

const GenericModel = Mimosa.AbstractMotifModel

function read_model(
    path::AbstractString, model_type::AbstractString="auto"; index::Integer=0, kwargs...
)
    format = model_type == "auto" ? :auto : Symbol(lowercase(model_type))
    return Mimosa.readmodel(path; format=format, index=index, kwargs...)
end

function read_model(path::AbstractString, model_type::Symbol; kwargs...)
    return read_model(path, String(model_type); kwargs...)
end

function write_model(path::AbstractString, model::Mimosa.AbstractMotifModel)
    return Mimosa.writemodel(path, model)
end

function rename_model(model::Mimosa.PWM, name::AbstractString)
    return Mimosa.PWM(name, model.representation, model.background)
end

function rename_model(model::Mimosa.BaMM, name::AbstractString)
    return Mimosa.BaMM(name, model.representation, model.order, model.motif_length)
end

function rename_model(model::Mimosa.SiteGA, name::AbstractString)
    return Mimosa.SiteGA(name, model.representation, model.motif_length)
end

function rename_model(model::Mimosa.Dimont, name::AbstractString)
    return Mimosa.Dimont(name, model.representation, model.span, model.motif_length)
end

function rename_model(model::Mimosa.Slim, name::AbstractString)
    return Mimosa.Slim(name, model.representation, model.span, model.motif_length)
end

motif_name(model::Mimosa.AbstractMotifModel) = String(Mimosa.modelname(model))
motif_length(model::Mimosa.AbstractMotifModel) = Mimosa.motif_length(model)
