$(document).ready(function () {
    // Function to load announcements and insert them into the DOM
    function loadAnnouncements(url, targetSelector) {
        $.get(url, function (data) {
            // Create a new div and add the fetched content
            const newDiv = $("<div>").html(data);

            // Insert the new div before the target element
            $(targetSelector).before(newDiv);
        }).fail(function () {
            console.error("Failed to load the message from the external file: " + url);
        });
    }

    // Load Immcantation announcements and insert into a new div
    loadAnnouncements("https://raw.githubusercontent.com/immcantation/immcantation/refs/heads/master/announcements.html", "#change-o-repertoire-clonal-assignment-toolkit");

    // Load this package announcements and insert into a new div
    loadAnnouncements("https://raw.githubusercontent.com/immcantation/changeo/refs/heads/master/announcements.html", "#change-o-repertoire-clonal-assignment-toolkit");    
});