#' Launch RShiny app for PPS analysis with mediation
#'
#' Launches an RShiny app to for PPS analysis. The
#' app is essentially a wrapper for the \code{pps}
#' function, with some additional visualizations.
#'
#'
#'
#'
#' @export


run.PPS.App.Mediation <- function() {

  ui <- shiny::fluidPage(

    # Application title
    shiny::titlePanel("PPS Analysis"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput("file", "Upload File", accept = ".csv"),

        shiny::uiOutput("penalty_choice"),
        shiny::uiOutput("node_1_choice"),
        shiny::uiOutput("node_2_choice"),
        shiny::uiOutput("mediator_choice"),
        shiny::uiOutput("K_choice"),
        shiny::uiOutput("button")


      ),


      shiny::mainPanel(
        shiny::tabsetPanel(type = "tabs",
                           shiny::tabPanel("Network", shiny::plotOutput("network_plot")),
                           shiny::tabPanel("Paths", shiny::verbatimTextOutput("paths")),
                           shiny::tabPanel("Mediation", shiny::verbatimTextOutput("mediation")),
                           shiny::tabPanel("Subnetwork", shiny::plotOutput("subnetwork")),
                           shiny::tabPanel("View Data", DT::dataTableOutput("dataDisplay"))

        )
      )
    )

  )

  # Define server logic required to draw a histogram
  server <- function(input, output) {


    output$penalty_choice <- shiny::renderUI({
      if (is.null(input$file$datapath)) {
        HTML("")
      } else {
        shiny::numericInput("penalty", "Graphical Lasso Penalty:", value = 1)
      }
    })

    output$node_1_choice <- shiny::renderUI({
      if (is.null(input$file$datapath)) {
        HTML("")
      }
      else {
        data <- read.csv(input$file$datapath)
        data <- data[,-1]
        names <- colnames(data)
        if (is.null(names)) {
          shiny::selectInput("node1", "Node 1", 1:dim(data)[2])
        } else {
          shiny::selectInput("node1", "Node 1", sort(names))
        }
      }
    })

    output$node_2_choice <- shiny::renderUI({
      if (is.null(input$file$datapath)) {
        HTML("")
      }
      else {
        data <- read.csv(input$file$datapath)
        data <- data[,-1]
        names <- colnames(data)
        if (is.null(names)) {
          shiny::selectInput("node2", "Node 2", 1:dim(data)[2])
        } else {
          shiny::selectInput("node2", "Node 2", sort(names))
        }
      }
    })


    output$mediator_choice <- shiny::renderUI({
      if (is.null(input$file$datapath)) {
        HTML("")
      }
      else {
        data <- read.csv(input$file$datapath)
        data <- data[,-1]
        names <- colnames(data)
        if (is.null(names)) {
          shiny::selectInput("mediator", "Mediator", 1:dim(data)[2])
        } else {
          shiny::selectInput("mediator", "Mediator", sort(names))
        }
      }
    })



    output$K_choice <- shiny::renderUI({
      if (is.null(input$file$datapath)) {
        HTML("")
      } else {
        shiny::numericInput("K",
                            "Search paths up to length: ",
                            value = 5,
                            min = 1,
                            max = NA,
                            step = 1)
      }
    })

    output$button <- shiny::renderUI({
      if (is.null(input$file$datapath)) {
        HTML("")
      } else {
        shiny::actionButton("button_submit", "Get PPS")
      }
    })

    #make everything update when submit is pressed
    file <- shiny::eventReactive(input$button_submit, {
      input$file
    })

    penalty <- shiny::eventReactive(input$button_submit, {
      input$penalty
    })

    node1 <- shiny::eventReactive(input$button_submit, {
      input$node1
    })

    node2 <- shiny::eventReactive(input$button_submit, {
      input$node2
    })

    mediator <- shiny::eventReactive(input$button_submit, {
      input$mediator
    })

    K <- shiny::eventReactive(input$button_submit, {
      input$K
    })


    output$network_plot <- shiny::renderPlot({
      if (is.null(input$file$datapath)) {
        plot.new()
      } else {

        data <- read.csv(input$file$datapath)
        data <- data[,-1]
        nnodes <- dim(data)[2]
        if (is.null(colnames(data))) {
          i <- input$node1
          j <- input$node2
          med <- input$mediator
        } else {
          i <- which(colnames(data) == input$node1)
          j <- which(colnames(data) == input$node2)
          med <- which(colnames(data) == input$mediator)
        }

        cov <- cov(as.matrix(unname(data)))


        pal <- c(igraph::categorical_pal(1), "blue", "red")

        res_gl <- glasso::glasso(cov, input$penalty)
        e <- ifelse(abs(res_gl$wi) > 0, 1, 0)
        diag(e) <- 0
        gr_gl <- igraph::graph_from_adjacency_matrix(e, mode = "undirected")

        temp <- rep(1, nnodes)
        temp[i] <- 2
        temp[j] <- 2
        temp[med] <- 3
        temp_size <- rep(5, nnodes)
        #temp_size[i] <- 7
        #temp_size[j] <- 7
        set.seed(1)

        plot(gr_gl, vertex.label.cex = 0.9,
             vertex.size = temp_size,
             vertex.label = names(data),
             vertex.color = pal[temp])
      }
    },
    height = 1000,
    width = 1000)

    output$paths <- shiny::renderPrint({
      data <- read.csv(file()$datapath)
      data <- data[,-1]
      if (floor(K()) != K()) {
        stop("Path length is not an integer.")
      }
      nnodes <- dim(data)[2]
      if (is.null(colnames(data))) {
        i <- node1()
        j <- node2()
      } else {
        i <- which(colnames(data) == node1())
        j <- which(colnames(data) == node2())
      }

      cov <- cov(as.matrix(unname(data)))



      res_gl <- glasso::glasso(cov, penalty())


      #convert precision matrix to partial correlation matrix
      pcor <- flip_off_diag(cov2cor(res_gl$wi))

      #add back names
      colnames(pcor) <- colnames(data)

      res <- pps(pcor, i, j, K = K(), prec = FALSE)



      print(Map(c, res$path, res$pps))


    })


    output$mediation <- shiny::renderPrint({
      data <- read.csv(file()$datapath)
      data <- data[,-1]
      if (floor(K()) != K()) {
        stop("Path length is not an integer.")
      }
      nnodes <- dim(data)[2]
      if (is.null(colnames(data))) {
        i <- node1()
        j <- node2()
        med <- mediator()
      } else {
        i <- which(colnames(data) == node1())
        j <- which(colnames(data) == node2())
        med <- which(colnames(data) == mediator())
        med <- colnames(data)[med]
      }



      cov <- cov(as.matrix(unname(data)))



      res_gl <- glasso::glasso(cov, penalty())


      #convert precision matrix to partial correlation matrix
      pcor <- flip_off_diag(cov2cor(res_gl$wi))

      #add back names
      colnames(pcor) <- colnames(data)

      res <- pps(pcor, i, j, K = K(), prec = FALSE)

      npaths <- length(res$path)

      pm <- 0
      for (h in 1:npaths) {
        if (med %in% res$path[[h]]) {
          pm <- pm + unlist(res$pps[[h]])
        }
      }

      print(pm)
      #get paths with mediator in them
      #f <-function(x) {
      #  med %in% x
      #}

      #inds <- lapply(res$path, f) %>% unlist()

      #print(sum(unlist(res$pps[inds])))

    })



    output$subnetwork <- shiny::renderPlot({
      data <- read.csv(file()$datapath)
      data <- data[,-1]
      if (floor(K()) != K()) {
        stop("Path length is not an integer.")
      }
      nnodes <- dim(data)[2]
      if (is.null(colnames(data))) {
        i <- node1()
        j <- node2()
      } else {
        i <- which(colnames(data) == node1())
        j <- which(colnames(data) == node2())
      }

      cov <- cov(as.matrix(unname(data)))



      res_gl <- glasso::glasso(cov, penalty())


      e <- ifelse(abs(res_gl$wi) > 0, 1, 0)
      diag(e) <- 0
      gr_gl <- igraph::graph_from_adjacency_matrix(e, mode = "undirected")
      V(gr_gl)$name <- colnames(data)


      #convert precision matrix to partial correlation matrix
      pcor <- flip_off_diag(cov2cor(res_gl$wi))

      #add back names
      colnames(pcor) <- colnames(data)

      res <- pps(pcor, i, j, K = K(), prec = FALSE, use.names = FALSE)

      #get indices appearing in top 10 paths

      a <- unique(unlist(res$path))

      #get subnetwork with only those indices



      #get edge ids for top path

      top_path <- res$path[[1]]
      ind_top_path <- top_path
      length_top_path <- length(top_path)
      if (length_top_path >= 3) {
        ind_top_path <- rep(0, 2*(length_top_path-2) + 2)
        for (l in 1:(length_top_path-1)) {
          ind_top_path[c((2*l-1),(2*l))] <- c(top_path[l], top_path[l+1])
        }
      }
      edges <- get.edge.ids(gr_gl, ind_top_path)


      set.seed(1)

      E(gr_gl)$color <- "gray"
      E(gr_gl)$width <- 1

      E(gr_gl)$width[edges] <- 2
      E(gr_gl)$color[edges] <- "black"

      pal <- c(igraph::categorical_pal(1), "blue")

      V(gr_gl)$color <- rep(pal[1], nnodes)
      V(gr_gl)$color[c(i,j)] <- rep(pal[2], 2)



      sub <- induced_subgraph(gr_gl, a)
      plot(sub,
           vertex.size = 8,
           #vertex.color = color_ac,
           #vertex.label = label_ac,
           vertex.label.color = "black",
           vertex.label.size = 4,
           vertex.label.family = "sans")

    },
    height = 1000,
    width = 1000)




    output$dataDisplay <- DT::renderDataTable({
      if (is.null(input$file$datapath)) {
        stop("No file detected yet.")
      }
      data <- read.csv(input$file$datapath)
      data <- data[,-1]
      data
    })



  }
  shiny::shinyApp(ui, server)
}
