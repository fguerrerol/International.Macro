

function plot_cons(sw::SOE; indiv=false)
jA, jz, Na = 5, 5,length(sw.agrid)
cons_mat = [sw.v[:,:,:][ja, jA, jz] for ja in eachindex(sw.agrid),
            jA in eachindex(sw.agrid)]
cons_agg = [sw.v[:,:,:][ja, ja, jz] for ja in eachindex(sw.agrid)]
colvec = [get(ColorSchemes.davos, (jA-1)/(Na-1)) for jA in eachindex(sw.agrid)]
scats = [scatter(x=sw.agrid, y=cons_mat[:, jA], 
        marker_color=colvec[jA], name ="A=$(@sprintf("%0.3g",Av))") for (jA, Av) in enumerate(sw.agrid)]


indiv || push!(scats, scatter(x=sw.agrid, y=cons_agg,
 line_dash="dash", line_width=3,name="Agregado", line_color="#710627"))

layout = Layout(title="Consumo",font_family ="Lato", 
        font_size = 18, width = 1920*0.5, height=1080*0.5,
        paper_bgcolor="#1e1e1e", plot_bgcolor="#1e1e1e",
         font_color="white",xaxis = attr(zeroline =false,
          gridcolor="#353535", title="<i>a"),
          yaxis = attr(zeroline =false, gridcolor="#353535"),
          )

plot(scats, layout)
end



