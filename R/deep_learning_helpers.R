plot_grayscale = function(image_array, label = NULL, class_names = NULL){
  if(!is.null(label) & !is.null(class_names)){
    title = class_names |> dplyr::filter(class == !!label) |> dplyr::pull(name)
  } else{
    title = label
  }
  p = image_array |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "row") |>
    dplyr::mutate_all(as.numeric) |>
    tidyr::pivot_longer(-row,
                 values_to = "intensity",
                 names_to = "col",
                 names_prefix = "V",
                 names_transform = list(col = as.integer)) |>
    dplyr::mutate(row = as.integer(row)) |>
    ggplot2::ggplot(ggplot2::aes(x = col, y = - row, fill = intensity)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(expand = c(0,0)) +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_fill_gradient(low = "white", high = "black") +
    ggplot2::theme_void() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.position = "none")

  if(!is.null(title)){
    p = p + ggplot2::ggtitle(title)
  }

  p
}

plot_confusion_matrix = function(predicted_responses, actual_responses, class_names = NULL){
  if(is.null(class_names)){
    rotate_angle = 0
    h_just = 0.5
    classes = sort(unique(c(predicted_responses, actual_responses)))
    class_names = tibble::tibble(class = classes, name = classes)
  } else{
    rotate_angle = 45
    h_just = 0
  }
  tibble::tibble(predicted_response = predicted_responses,
         actual_response = actual_responses) |>
    table() |>
    tibble::as_tibble() |>
    ggplot2::ggplot(ggplot2::aes(x = predicted_response, y = actual_response, label = n, fill = n)) +
    ggplot2::geom_tile(alpha = 0.5) +
    ggplot2::geom_text() +
    ggplot2::scale_x_discrete(position = "top",
                     expand = c(0,0),
                     breaks = class_names$class,
                     labels = class_names$name) +
    ggplot2::scale_y_discrete(limits = rev,
                     expand = c(0,0),
                     breaks = class_names$class,
                     labels = class_names$name) +
    ggplot2::scale_fill_gradient(low = "blue", high = "red") +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "Predicted Response",
         y = "Actual Response") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = rotate_angle, vjust = 0.5, hjust=h_just))
}

plot_model_history = function(history){
  epochs = 1:length(history$loss)
  dplyr::bind_rows(tibble::tibble(epoch = epochs,
                   metric = "loss",
                   value = history$loss,
                   set = "training"),
            tibble::tibble(epoch = epochs,
                   metric = "accuracy",
                   value = history$accuracy,
                   set = "training"),
            tibble::tibble(epoch = epochs,
                   metric = "loss",
                   value = history$val_loss,
                   set = "validation"),
            tibble::tibble(epoch = epochs,
                   metric = "accuracy",
                   value = history$val_accuracy,
                   set = "validation")
  ) |>
    dplyr::mutate(metric = factor(metric, levels = c("loss", "accuracy"))) |>
    ggplot2::ggplot(ggplot2::aes(x = epoch, y = value, colour = set)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_grid(metric ~ ., scales = "free") +
    ggplot2::labs(x = "Epoch") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
          legend.title = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank())
}
