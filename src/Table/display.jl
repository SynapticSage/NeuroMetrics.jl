
function display_table(data; others...)
    w=Window()
    body!(w, TableView.showtable(data, dark=true, height=950, others...))
end

