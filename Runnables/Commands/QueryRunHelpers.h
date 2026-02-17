#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include "../../DataStructures/RAPTOR/Data.h"

template<typename Algorithm, typename Queries>
void runRAPTORBasedQueriesAndWriteJourneyStats(Algorithm &algorithm, const Queries &queries,
                                               const std::string &outFileName,
                                               const std::string &detailedOutFileName) noexcept {
    const size_t n = queries.size();
    std::cout << "Running queries ..." << std::flush;
    ProgressBar progressBar(n);
    progressBar.SetDotOutputStep(1);
    progressBar.SetPercentOutputStep(5);
    double numJourneys = 0;

    std::ofstream overviewOut(outFileName);
    if (!overviewOut.good()) {
        std::cerr << "Could not open output file " << outFileName << " for writing journeys.";
        return;
    }
    overviewOut << "request_id,journey,departure_time,arrival_time,num_trips,"
            "accegr_transfer_time,intermediate_transfer_time,wait_time,in_vehicle_time\n";

    std::ofstream detailedOut(detailedOutFileName);
    if (!detailedOut.good()) {
        std::cerr << "Could not open output file " << detailedOutFileName << " for writing detailed journeys.";
        return;
    }
    detailedOut << "request_id,journey,leg,from,to,departure_time,arrival_time,type\n";

    size_t requestId = 0;
    for (const auto &query: queries) {
        algorithm.run(query.source, query.departureTime, query.target);
        numJourneys += algorithm.getJourneys().size();
        const auto journeys = algorithm.getJourneys();
        const auto arrivals = algorithm.getArrivals();
        Assert(journeys.size() == arrivals.size(), "Number of journeys and arrivals do not match.");
        for (size_t i = 0; i < journeys.size(); ++i) {
            const RAPTOR::Journey &journey = journeys[i];
            const auto &arrival = arrivals[i];
            const int accEgrTransferTime = RAPTOR::initialTransferTime(journey);
            const int intermediateTransferTime = RAPTOR::intermediateTransferTime(journey);
            int inVehicleTime = 0;
            for (const auto &leg: journey) {
                if (leg.usesRoute) {
                    inVehicleTime += leg.arrivalTime - leg.departureTime;
                }
            }
            const int waitTime = arrival.arrivalTime - query.departureTime - inVehicleTime - accEgrTransferTime -
                                 intermediateTransferTime;
            const auto journeyId = i;
            overviewOut << requestId << "," << journeyId << "," << query.departureTime << ","
                    << arrival.arrivalTime << "," << arrival.numberOfTrips << "," << accEgrTransferTime << ","
                    << intermediateTransferTime << "," << waitTime << "," << inVehicleTime << "\n";

            for (size_t j = 0; j < journey.size(); ++j) {
                const auto &leg = journey[j];
                const auto legId = j;
                detailedOut << requestId << "," << journeyId << "," << legId << "," << leg.from.value() << ","
                        << leg.to.value() << "," << leg.departureTime << "," << leg.arrivalTime << ","
                        << (leg.usesRoute ? "route" : "transfer") << "\n";
            }
        }
        ++requestId;
        ++progressBar;
    }
    overviewOut.close();
    detailedOut.close();
    std::cout << " done." << std::endl;
    algorithm.getProfiler().printStatistics();
    std::cout << "Avg. journeys: " << String::prettyDouble(numJourneys / n) << std::endl;
}

template<typename Algorithm, typename Queries>
void runMcRAPTORBasedQueriesAndWriteJourneyStats(Algorithm &algorithm, const Queries &queries,
                                                 const std::string &outFileName,
                                                 const std::string &detailedOutFileName) noexcept {
    const size_t n = queries.size();
    std::cout << "Running queries ..." << std::flush;
    ProgressBar progressBar(n);
    progressBar.SetDotOutputStep(1);
    progressBar.SetPercentOutputStep(5);
    double numJourneys = 0;

    std::ofstream overviewOut(outFileName);
    if (!overviewOut.good()) {
        std::cerr << "Could not open output file " << outFileName << " for writing journeys.";
        return;
    }
    overviewOut << "request_id,journey,departure_time,arrival_time,num_trips,"
            "accegr_transfer_time,intermediate_transfer_time,wait_time,in_vehicle_time\n";

    std::ofstream detailedOut(detailedOutFileName);
    if (!detailedOut.good()) {
        std::cerr << "Could not open output file " << detailedOutFileName << " for writing detailed journeys.";
        return;
    }
    detailedOut << "request_id,journey,leg,from,to,departure_time,arrival_time,type\n";

    size_t requestId = 0;
    for (const StopQuery &query: queries) {
        algorithm.run(query.source, query.departureTime, query.target);
        numJourneys += algorithm.getJourneys().size();
        const auto journeys = algorithm.getJourneys();
        const auto paretoLabels = algorithm.getResults();
        Assert(journeys.size() == paretoLabels.size(), "Number of journeys and paretoLabels do not match.");
        for (size_t i = 0; i < journeys.size(); ++i) {
            const RAPTOR::Journey &journey = journeys[i];
            const auto &label = paretoLabels[i];
            const int accEgrTransferTime = RAPTOR::initialTransferTime(journey);
            const int intermediateTransferTime = RAPTOR::intermediateTransferTime(journey);
            Assert(accEgrTransferTime + intermediateTransferTime == label.walkingDistance,
                   "Transfer times from journey do not match walking distance from label.");
            int inVehicleTime = 0;
            for (const auto &leg: journey) {
                if (leg.usesRoute) {
                    inVehicleTime += leg.arrivalTime - leg.departureTime;
                }
            }
            const int waitTime = label.arrivalTime - query.departureTime - inVehicleTime - accEgrTransferTime -
                                 intermediateTransferTime;
            const auto journeyId = i;
            overviewOut << requestId << "," << journeyId << "," << query.departureTime << ","
                    << label.arrivalTime << "," << label.numberOfTrips << "," << accEgrTransferTime << ","
                    << intermediateTransferTime << "," << waitTime << "," << inVehicleTime << "\n";

            for (size_t j = 0; j < journey.size(); ++j) {
                const auto &leg = journey[j];
                const auto legId = j;
                detailedOut << requestId << "," << journeyId << "," << legId << "," << leg.from.value() << ","
                        << leg.to.value() << "," << leg.departureTime << "," << leg.arrivalTime << ","
                        << (leg.usesRoute ? "route" : "transfer") << "\n";
            }
        }
        ++requestId;
        ++progressBar;
    }
    overviewOut.close();
    detailedOut.close();
    std::cout << " done." << std::endl;
    algorithm.getProfiler().printStatistics();
    std::cout << "Avg. journeys: " << String::prettyDouble(numJourneys / n) << std::endl;
}
