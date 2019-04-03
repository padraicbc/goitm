package main

import (
	"flag"
	"fmt"

	"../../pkg/convert"
)

func main() {
	var xe, yn float64

	flag.Float64Var(&xe, "x", 0, "east")
	flag.Float64Var(&yn, "y", 0, "north")

	flag.Parse()
	fmt.Println(convert.Itm2wgs84(yn, xe))

}
