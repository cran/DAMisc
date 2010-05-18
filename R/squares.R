squares <- function(ll, width=1,col){ 
    poly.x <- c(ll[1], ll[1]+width, ll[1]+width, ll[1], ll[1])
    poly.y <- c(ll[2], ll[2], ll[2]+width, ll[2]+width, ll[2])
    polygon(poly.x, poly.y, col=col)
}