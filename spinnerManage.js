$(document).ready(function() {
            $('#klickButton').click(function() {
            $(".loading-spinner").show();
        });  
    });
    $(document).on("shiny:connected", function(e) {
            $(".loading-spinner").hide();
    });
